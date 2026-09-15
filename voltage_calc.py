#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 13:46:54 2020
new version in 1.4 branch 19 May 2024. No substantive changes to this file

@author: perkel
"""

# implant voltage calculations from Goldwyn et al, 2010, denovo.

# Import dependencies
import cProfile
import io
import os
import pickle
import pstats
import time
import concurrent.futures as cf
from datetime import datetime
from pstats import SortKey
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integ
import scipy.special as spec
import remove_nans as r_nans


def goldwyn_beta(eps, k2, rs, re, n):
    # 'n' may be a scalar or an array of orders -- scipy.special's Bessel functions broadcast
    # over the order argument, so passing the full array of orders at once (see integ_func below)
    # computes all of them in one vectorized call instead of one Python-level call per order.
    #
    # Numerical note: for small k2 (the quadrature variable can get very close to the lower
    # integration bound, ~1e-12) combined with high order n, I_n(krs) underflows to exactly 0.0
    # while K_n(krs)/K_n'(krs) overflow to +-inf (small-argument Bessel asymptotics: I_n(z) ~
    # (z/2)^n/n! -> 0, K_n(z) ~ (n-1)!/2 * (2/z)^n -> inf as z -> 0, and both blow up/vanish faster
    # for larger n). That makes 'denom' below hit an IEEE inf*0 = nan for the highest orders in
    # exactly this regime -- confirmed to happen in real runs (not just a theoretical edge case),
    # e.g. order 46 at krs ~ 8.3e-6. Where I_n(krs) has genuinely underflowed to 0.0, the true
    # physical contribution of that order is astronomically small, so we force it to exactly 0
    # instead of letting the inf*0 corruption through -- verified to reproduce the correct,
    # smoothly-varying result (checked against neighboring z values) with zero effect anywhere
    # the computation was already stable.
    krs = k2 * rs
    with np.errstate(invalid='ignore'):
        iv_krs = spec.iv(n, krs)
        denom = ((eps * spec.kvp(n, krs) * iv_krs) - (spec.kn(n, krs) * spec.ivp(n, krs)))
        result = spec.iv(n, k2 * re) / (denom * krs)
    negligible = (iv_krs == 0.0)
    return result, negligible


def goldwyn_phi(eps, k, rs, re, reval, n):  # Note Goldwyn alpha is zero for eval pt. outide cylinder
    beta, negligible = goldwyn_beta(eps, k, rs, re, n)
    with np.errstate(invalid='ignore'):
        phi = beta * spec.kn(n, k * reval)
    # See the note in goldwyn_beta: force genuinely-negligible orders to exactly 0 rather than
    # propagating a NaN from either this multiplication or the one inside goldwyn_beta.
    return np.where(negligible, 0.0, phi)


def integ_func(x, m_max, pratio, rad, reval, z, theta, relec):  # This is the Bessel function along z axis
    # Vectorized over the order index (0..m_max-1): previously this ran m_max separate Python-level
    # calls to goldwyn_phi/goldwyn_beta (one per order), each doing its own scipy.special dispatch.
    # Since those functions broadcast over an array-valued order, computing all m_max orders in one
    # call is equivalent (verified numerically identical to the old per-order loop, to within
    # floating-point noise) and substantially faster -- this function is called many times per
    # integ.quad() evaluation, across thousands of (electrode, z, y) combinations.
    n_arr = np.arange(m_max)
    gamma = np.full(m_max, 2.0)
    gamma[0] = 1.0
    phi_vals = goldwyn_phi(pratio, x, rad, relec, reval, n_arr)
    sum_contents = np.sum(gamma * np.cos(n_arr * theta) * phi_vals)
    # (The old code also tracked a per-order 'rel_incrs' relative-increment array for debugging;
    # it was never read after being computed, so it's dropped here rather than vectorized too.)

    return np.cos(x * z) * sum_contents


def _process_electrode(i, r_elec, z_eval, m_max, res_ratio, cyl_radius, reval_base, int_start, int_end,
                        itol, y_vals, n_yeval, n_ypos, y_inc, n_relec):
    """Worker task: compute the full (y, z) voltage/activation grid for one electrode radial
    position. Electrode positions are fully independent of each other (each only reads its own
    'r_elec' and writes its own row of voltage_vals/activation_vals), so this loop -- previously
    39 serial iterations -- is a natural fit for a process pool."""
    n_z = len(z_eval)
    voltage_row = np.zeros((n_ypos, n_z))
    activation_row = np.zeros(n_z)
    trans_voltage = np.empty(n_ypos)  # private to this worker; the original shared one array across electrodes

    print('starting electrode position ', i, ' of ', n_relec)

    for m, this_z in enumerate(z_eval):
        # loop on y positions to get 2nd spatial derivative
        for j, y_val in enumerate(y_vals[0:n_yeval + 1]):
            this_theta = np.arctan(y_val / reval_base)
            r_prime = np.sqrt((np.power(y_val, 2)) + (reval_base ** 2))  # distance of eval pt. from center of cylinder
            [itemp, error] = integ.quad(integ_func, int_start, int_end, epsabs=itol, limit=1000,
                                        args=(m_max, res_ratio, cyl_radius, r_prime, this_z, this_theta, r_elec))
            temp_v = itemp / (2 * (np.pi ** 2))  # From Goldwyn eqn. 11
            voltage_row[j, m] = temp_v
            voltage_row[n_ypos - (j + 1), m] = temp_v  # place same value in mirror-symmetric position
            # Extra variable to make derivative calculation easy
            trans_voltage[j] = temp_v
            trans_voltage[n_ypos - (j + 1)] = temp_v

        # Calculate the second spatial derivative in the y dimension
        trans_v_prime = np.diff(np.diff(trans_voltage)) / (y_inc ** 2)
        activation_row[m] = trans_v_prime[n_yeval - 1]  # Value in center
        if np.isnan(activation_row[m]):
            print('nan value for i == ', i, " and m == ", m)

    return i, voltage_row, activation_row


if __name__ == '__main__':
    # See Rubinstein dissertation, 1988, Ch 6.

    run_start_time = time.perf_counter()  # wall-clock timer for the entire run (see the end of this
    # block for where it's reported). This is a proper wall-clock measurement covering the whole run,
    # including the parallelized electrode-position loop -- unlike the cProfile-based 'run_duration'
    # previously stored below, which only measures CPU time in the main process and therefore massively
    # under-counts the real time once the heavy computation moved into worker processes.

    # Main parameters to vary
    radius = 1.0  # cylinder radius
    res_int = 70.0  # internal resistivity
    res_ext = 250.0  # external resistivity
    vtable_dir = 'v_tables/'
    if not os.path.exists(vtable_dir):
        os.makedirs(vtable_dir)

    date = '14September2026'
    resol = '_MedResolution_'
    rext_text = 'Rext' + str(round(res_ext))
    r_int_text = '_Rint' + str(round(res_int)) + '.dat'
    output_filename = vtable_dir + date + resol + rext_text + r_int_text

    pr = cProfile.Profile()

    # Field parameters. zEval can be higher (more precise, slower) or lower resolution (less precise, faster)
    fp = {'model': 'cylinder_3d', 'fieldtype': 'Voltage_Activation', 'evaltype': 'SG', 'cylRadius': radius,
          'relec': np.arange(-0.95, 0.951, 0.05), 'resInt': res_int, 'resExt': res_ext, 'rspace': 1.3, 'theta': 0.0,
          # "High" resolution
          # 'zEval': (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,
          #              1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
          #              3.0, 3.2, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7,
          #              4.8, 4.9, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8,
          #              8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
          #              13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0,
          #              22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0),
          # Medium resolution
          'zEval': (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
                   8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
                   32.0, 34.0, 36.0, 38.0, 40.0),
          # Low resolution
          # 'zEval': (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
          #            1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
          #           8.5, 9.0, 9.5, 10.0, 15, 20.0, 25.0, 30.0, 35.0, 40),
          # Very low resolution for debugging
          #'zEval': (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
          'mMax': 47, 'intStart': 1e-12, 'intEnd': 500.0, 'reval': 1.3,
          'ITOL': 1e-6, 'runDate': 'rundate', 'runOutFile': 'savefile', 'run_duration': 0.0}

    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    fp['runDate'] = date_time
    fp['runOutFile'] = output_filename

    # Derived values
    res_ratio = fp['resInt'] / fp['resExt']
    n_z = len(fp['zEval'])  # this may be a problem if zEval is just a scalar
    n_ypos = 3  # # of values to calculate across the y dimension for the 2nd spatial derivative

    relec_range = fp['relec']  # electrode radial positions
    n_relec = len(relec_range)
    voltage_vals = np. zeros((n_relec, n_ypos, n_z))
    activation_vals = np.zeros((n_relec, n_z))
    if_plot = False

    pr.enable()  # Start the profiler
    # Note: this only profiles the main process below (task submission/collection), not the actual
    # per-electrode computation, which now runs in worker processes and won't show up here.

    n_yeval = n_ypos // 2  # floor division
    y_inc = 0.001  # 10 microns
    y_vals = np.arange(-n_yeval * y_inc, n_yeval * (y_inc * 1.01), y_inc)

    # loop on electrode radial positions -- each position is fully independent (reads only its own
    # r_elec, writes only its own row of voltage_vals/activation_vals), so this is parallelized
    # across processes rather than run serially.
    n_workers = max(1, (os.cpu_count() or 2) - 1)
    print('Parallelizing electrode-position loop across ', n_workers, ' worker processes')

    with cf.ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(_process_electrode, i, r_elec, fp['zEval'], fp['mMax'], res_ratio,
                                    fp['cylRadius'], fp['reval'], fp['intStart'], fp['intEnd'], fp['ITOL'],
                                    y_vals, n_yeval, n_ypos, y_inc, n_relec)
                   for i, r_elec in enumerate(relec_range)]

        n_done = 0
        for fut in cf.as_completed(futures):
            i, voltage_row, activation_row = fut.result()
            voltage_vals[i, :, :] = voltage_row
            activation_vals[i, :] = activation_row
            n_done += 1
            print('electrode position ', n_done, ' of ', n_relec, ' complete (i=', i, ')')

    pr.disable()  # stop the profiler

    # Display profiler results
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(10)
    print(s.getvalue())

    # Save the data into a single file
    # Note: ps.total_tt (cProfile's cumulative CPU time) only covers the main process -- the actual
    # per-electrode computation above now runs in worker processes and isn't included in it, so it
    # would badly under-report the real time since parallelization. Use the wall-clock timer instead.
    fp['run_duration'] = time.perf_counter() - run_start_time  # Place run duration in field parameters
    with open(fp['runOutFile'], 'wb') as combined_data:
        pickle.dump([fp, voltage_vals, activation_vals], combined_data)  # Save all data
    combined_data.close()

    nan_locs = np.argwhere(np.isnan(activation_vals))
    if nan_locs.size != 0:
        print('WARNING: There were NaN values in the act_vals array at locations: ', nan_locs)
        r_nans.remove_nans(fp['runOutFile'])
        print('These have been replaced by interpolated values in a new file with no_nans appended to the filename.')
    else:
        print('Complete: No nans are in the file.')

    run_elapsed_time = time.perf_counter() - run_start_time
    print(f'Total run duration: {run_elapsed_time:.1f} s ({run_elapsed_time / 60:.2f} min)')

    if if_plot:
        # Some plots that might be helpful
        #  plt.plot(relec_range, voltageVals[:, 10], 'or')
        #  plt.plot(relec_range, activationVals[:, 0], 'ob')
        #  plt.xlabel('Electrode radial position (mm)')
        #  plt.ylabel('Red: voltage; Blue: activation')
        #  plt.xlim(-1.0, 1.0)
        #  plt.yscale('log')

        plt.figure()
        plt.xlabel('y position (mm)')
        plt.ylabel('Voltage')
        plt.plot(y_vals, voltage_vals[0, :, 0], '.b')
        plt.plot(y_vals, voltage_vals[1, :, 0], '.g')
        plt.show()
