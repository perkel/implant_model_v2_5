# InverseModelCombined.py
# Script to fit TP and MP threshold data to the Goldwyn, Bierer, and
# Bierer cochlear activation model.

# Original code by Steve Bierer.

# Translated to python 3 and adapted to fit both survival and rpos values by David J. Perkel December 2020

# Modified 6 August 2021 by DJP to start by fitting a single electrode at a time. This gives initial conditions
# for each electrode. Then there's a holistic fitting process using all electrode parameters
# Latest version 1 December 2023.

import cProfile
import io
import pstats
from pstats import SortKey
import csv
import os
import scipy.signal as sig
import scipy.stats as stats
import matplotlib.pyplot as plt
import pickle
import scipy.optimize as opt
from lmfit import Minimizer, Parameters, report_fit
import intersection as intsec

# local files
import set_scenario as s_scen
import load_fwd_csv_data as lcsv
import surv_full
import get_thresholds as gt
import subject_data
import plot_inverse_results
from common_params import *  # import common values across all models

# User adjustable parameters
fit_mode = 'combined'  # Which variable(s) to fit? Alternatives are 'combined', 'rpos' or 'survival'
fit_using = 'mp_tp'  # Options are 'mp_tp' for thresholds, 'mp_mlevel' for mp threshuold and m-level,
# and 'ecap' for ecap growth curve threshold and slope
ifPlot = True  # Whether to plot output at end
unsupervised = True  # Makes & saves summary plots but does not display them and wait for user input before proceeding
ifPlotGuessContours = False  # Option to plot initial guesses for parameters given to the fitting algorithm
fit_tol = 0.1  # Fit tolerance for subject fits
use_minimizer = True  # If true, uses the lmfit wrapper around scipy optimize. Otherwise, vanilla scipy.minimize
if_save_npy = True
m_test = 'False'


# For optimizing fit to thresholds need e_field, sim_params, sigvals
# This is for all electrodes at once


def objectivefunc_lmfit_all(par, sigvals, sim_params, f_par, e_field, thr_goals):
    # Repack parameters into arrays
    nel = len(sim_params['electrodes']['zpos'])
    vals = par.valuesdict()
    show_retval = True  # helpful to track progress of the fitting process

    sim_params['electrodes']['rpos'] = np.zeros(nel)
    if not tp_extend:
        for i in range(0, nel - 2):
            varname = 'v_%i' % i
            myvalue = vals[varname]
            sim_params['electrodes']['rpos'][i + 1] = myvalue
            sim_params['electrodes']['rpos'][0] = sim_params['electrodes']['rpos'][1]
            sim_params['electrodes']['rpos'][-1] = sim_params['electrodes']['rpos'][-2]

            tempsurv = np.zeros(nel - 2)
            tempsurv2 = np.zeros(nel)
            for j, loopval in enumerate(range(nel - 2, 2 * (nel - 2))):
                varname = 'v_%i' % (j + nel - 2)
                # print('i: ', i, ' loopval: ', loopval, ' varname: ', varname)
                myvalue = vals[varname]
                tempsurv[j] = myvalue
            tempsurv2[1:nel - 1] = tempsurv
            tempsurv2[0] = tempsurv2[1]
            tempsurv2[-1] = tempsurv2[-2]
            sim_params['neurons']['nsurvival'] = surv_full.surv_full(sim_params['electrodes']['zpos'],
                                                                     tempsurv2, simParams['grid']['z'])

    else:
        for i in range(0, nel):
            varname = 'v_%i' % i
            myvalue = vals[varname]
            sim_params['electrodes']['rpos'][i] = myvalue

        tempsurv = np.zeros(nel)
        for i, loopval in enumerate(range(nel, 2 * nel)):
            varname = 'v_%i' % (i + nel)
            myvalue = vals[varname]
            tempsurv[i] = myvalue

            tempsurv[0] = tempsurv[1]
            tempsurv[-1] = tempsurv[-2]
        sim_params['neurons']['nsurvival'] = surv_full.surv_full(sim_params['electrodes']['zpos'],
                                                                 tempsurv, simParams['grid']['z'])

    # Call for monopolar then tripolar
    sim_params['channel']['sigma'] = sigvals[0]
    sim_params['channel']['behavior'] = 0
    thresh_mp = gt.get_thresholds(e_field, f_par, sim_params)

    if fit_using == 'mp_tp':
        sim_params['channel']['sigma'] = sigvals[1]   # 0 if looking at M levels
        sim_params['channel']['behavior'] = 0
        thresh_tp = gt.get_thresholds(e_field, f_par, sim_params)
    elif fit_using == 'mp_mlevel':
        sim_params['channel']['behavior']=1  # if looking at M levels
        thresh_tp = gt.get_thresholds(e_field, f_par, sim_params)
    elif fit_using == 'ecap':
        # not sure what goes here yet
        # TODO put ecap code here
        pass
    else:
        pass


    # Calculate errors
    mp_err = np.nanmean(np.abs(np.subtract(thresh_mp[0], thr_goals['thrmp_db'])))
    tp_err = np.nanmean(np.abs(np.subtract(thresh_tp[0][1:nel - 1], thr_goals['thrtp_db'][1:nel - 1])))
    mean_error = (mp_err + tp_err) / 2.0
    mp_diff = np.subtract(thresh_mp[0], thr_goals['thrmp_db'])
    tp_diff = np.subtract(thresh_tp[0][1:nel - 1], thr_goals['thrtp_db'][1:nel - 1])

    tempzero = np.zeros(1)
    if tp_extend:
        retval = np.concatenate((mp_diff, tempzero, tp_diff, tempzero))
    else:
        retval = np.concatenate((mp_diff[1:-1], tp_diff))
    # Returns a vector of errors, with the first and last of the tripolar errors set to zero
    # because they can't be calculated
    if show_retval:  # helpful for debugging
        scen = simParams['run_info']['scenario']
        print('subj/scen: ', scen, '; tempsurv[4] = ', '%.3f' % float(tempsurv[4]), '; rpos[4]= ',
              '%.3f' % sim_params['electrodes']['rpos'][4],
              '; Mean abs error (dB) = ', '%.3f' % mean_error, '; Max error (dB) = ',
              '%.3f' % np.nanmax(np.abs(retval)))
    return retval


# For optimizing fit to thresholds need e_field, sim_params, sigvals
# This is for all electrodes at once
def objectivefunc_minimize_all(par, sigvals, sim_params, f_par, e_field, thr_goals):
    # Repack parameters into arrays
    nel = len(sim_params['electrodes']['zpos'])
    show_retval = True  # helpful to track progress of the fitting process

    sim_params['electrodes']['rpos'] = par[0:nel]
    tempsurv = par[nel:]
    sim_params['neurons']['nsurvival'] = surv_full.surv_full(sim_params['electrodes']['zpos'],
                                                             tempsurv, simParams['grid']['z'])

    # Call for monopolar then tripolar
    sim_params['channel']['sigma'] = 0  # sigvals[0]
    sim_params['channel']['behavior'] = 0
    thresh_mp = gt.get_thresholds(e_field, f_par, sim_params)
    sim_params['channel']['sigma'] = sigvals[1]  # 0 for m level
    sim_params['channel']['behavior'] = 0  # 1 for m level
    # sim_params['channel']['behavior']=1 for target changing w m levels
    thresh_tp = gt.get_thresholds(e_field, f_par, sim_params)
    # Calculate errors
    mp_err = np.nanmean(np.abs(np.subtract(thresh_mp[0], thr_goals['thrmp_db'])))
    tp_err = np.nanmean(np.abs(np.subtract(thresh_tp[0][1:nel - 1], thr_goals['thrtp_db'][1:nel - 1])))
    mean_error = (mp_err + tp_err) / 2.0
    # retval = np.append(np.subtract(thresh_mp[0], thr_goals['thrmp_db']), 0)
    mp_diff = np.subtract(thresh_mp[0], thr_goals['thrmp_db'])
    tp_diff = np.subtract(thresh_tp[0][1:nel - 1], thr_goals['thrtp_db'][1:nel - 1])
    tempzero = np.zeros(1)
    retval = np.concatenate((mp_diff, tempzero, tp_diff, tempzero))
    # Returns mean error, with the first and last of the tripolar errors set to zero
    # because they can't be calculated
    if show_retval:  # helpful for debugging
        scen = simParams['run_info']['scenario']
        print('subj/scen: ', scen, '; tempsurv[4] = ', '%.3f' % tempsurv[4], '; rpos[4]= ',
              '%.3f' % sim_params['electrodes']['rpos'][4],
              '; Mean abs error (dB) = ', '%.3f' % mean_error, '; Max error (dB) = ',
              '%.3f' % np.nanmax(np.abs(retval)))
    return mean_error


def max_diff_adjacent(pars):
    # pars 0-15 are for distance
    diff = np.diff(pars[0:16])
    abs_diff = np.abs(diff)
    max_abs_diff = np.max(abs_diff)
    return max_abs_diff


def find_closest(x1, y1, x2, y2):  # returns indices of the point on each curve that are closest
    # Brute force (pretty ugly, but hopefully a rare case)
    # TODO Should test for case where x1 or x2 is a scalar
    n1 = len(x1)
    n2 = len(x2)

    min_dist = 5.0
    min_idx = [0, 0]
    for ii in range(n1):
        for jj in range(n2):
            dist = np.sqrt((((x1[ii] - x2[jj]) / 2) ** 2) + ((y1[ii] - y2[jj]) ** 2))
            if dist < min_dist:
                min_idx = [ii, jj]
                min_dist = dist

    return min_idx


def is_scenario(scen):  # test whether this is an artificial scenario or subject
    # if this scenario is a subject, set use_forward_model to be false
    if (scen[0] == 'A' or scen[0] == 'S') and scen[1:3].isnumeric():
        return False  # it's a subject
    else:
        return True


def inverse_model_combined(mode):  # Start this script

    if mode == 'main':
        pass  # proceed as normal
    elif mode == 'survey':
        # load survey params
        param_file = 'surv_params.txt'
        tempdata = np.zeros(4)  # 4 values
        if os.path.exists(param_file):
            with open(param_file, newline='') as csvfile:
                datareader = csv.reader(csvfile, delimiter=',')
                ncol = len(next(datareader))
                csvfile.seek(0)
                for i, row in enumerate(datareader):
                    # Do the parsing
                    tempdata[i] = row[0]

        res_ext = tempdata[0]
        NEURONS['act_stdrel'] = tempdata[1]
        # NEURONS['thrtarg'] = tempdata[2]
        NEURONS['act_ctr'] = tempdata[2]
        NEURONS['thrtarg'] = 100
        espace = tempdata[3]
    else:  # should not happen
        print('fwd_model called with unrecognized mode: ', mode)
        exit()

    if not os.path.isdir(INV_OUT_PRFIX):
        os.mkdir(INV_OUT_PRFIX)
        os.mkdir(INVOUTPUTDIR)
    else:
        if not os.path.isdir(INVOUTPUTDIR):
            os.mkdir(INVOUTPUTDIR)

    # First make sure that the 2D forward model has been run and load the data
    # It would be ideal to double check that it's the correct 2D data with identical parameters
    pr = cProfile.Profile()
    pr.enable()  # Start the profiler

    # Open field table and load data
    if "fieldTable" not in locals():
        with open(FIELDTABLE, 'rb') as combined_data:
            data = pickle.load(combined_data)
            combined_data.close()
            # (output is 'fieldTable' and 'fieldParams')
            fp = data[0]
            # Temp fixup
            fp['zEval'] = np.array(fp['zEval'])
            act_vals = data[2]  # the data[1] has voltage values, which we are ignoring here
            simParams['grid']['table'] = act_vals

    # Now hold these data until about to fit a particular combination of monopolar and tripolar threshold values
    # to see if there is more than one solution

    surv_grid_vals = np.arange(0.04, 0.97, 0.02)
    # TODO should read these from the file, not have the range hard coded here
    # TODO should we change these grid values to restrict the position a bit more??
    rpos_grid_vals = np.arange(-0.95, 0.96, 0.02)

    n_elec_pos = 0
    n_z_pos = 0

    num_scen = len(scenarios)
    # Set up array for summary values
    thresh_err_summary = np.zeros((num_scen, 2))
    thresh_summary_full = np.zeros((num_scen, 2, 2, NELEC))  # indices: scenario, (target, fit), sigma, electrodes
    rpos_summary = []
    rpos_err_summary = np.zeros(num_scen)
    surv_err_summary = np.zeros(num_scen)
    dist_corr = np.zeros(num_scen)
    dist_corr_p = np.zeros(num_scen)

    for scen in range(0, len(scenarios)):
        scenario = scenarios[scen]

        use_fwd_model = is_scenario(scenario)

        simParams['run_info']['scenario'] = scenario

        if use_fwd_model:
            [survVals, electrodes['rpos'], espace] = s_scen.set_scenario(scenario, NELEC)
        else:
            subject = scenario
            retval = subject_data.subj_thr_data(subject)
            thr_data = {'thrmp_db': retval[0], 'thrmp': [], 'thrtp_db': retval[1], 'thrtp': [], 'thrtp_sigma': 0.9}
            sigVals = retval[2]
            espace = retval[3]

        # Load 2D results corresponding to the correct electrode spacing
        if espace == 0.85:
            e_txt = '085'
        elif espace == 1.1:
            e_txt = '110'
        else:
            e_txt = 'xxx'
        es_text = '_espace_' + e_txt

        datafile = FWDOUTPUTDIR + "Monopolar_2D_" + STD_TEXT + es_text + ".csv"
        file = open(datafile)
        numlines = len(file.readlines())
        file.close()

        with open(datafile, newline='') as csvfile:
            datareader = csv.reader(csvfile, delimiter=',')
            ncol = len(next(datareader))
            csvfile.seek(0)
            mono_thr = np.empty([numlines, ncol])
            for i, row in enumerate(datareader):
                # Do the parsing
                mono_thr[i, :] = row

        # Load tripolar data
        datafile = FWDOUTPUTDIR + "Tripolar_09_2D_" + STD_TEXT + es_text + ".csv"
        file = open(datafile)
        numlines = len(file.readlines())
        file.close()

        with open(datafile, newline='') as csvfile:
            datareader = csv.reader(csvfile, delimiter=',')
            ncol = len(next(datareader))
            csvfile.seek(0)
            tripol_thr = np.empty([numlines, ncol])
            for i, row in enumerate(datareader):
                # Do the parsing
                tripol_thr[i, :] = row

        ELEC_MIDPOINT = GRID['z'][-1] / 2.0  # electrode array midpoint
        ARRAY_BASE = -(np.arange(NELEC - 1, -1, -1) * espace)
        array_mid = (ARRAY_BASE[0] + ARRAY_BASE[-1]) / 2.0
        electrodes['zpos'] = (ELEC_MIDPOINT - array_mid) + ARRAY_BASE

        if use_fwd_model:
            csv_file = FWDOUTPUTDIR + 'FwdModelOutput_' + scenario + '.csv'
            [thr_data, ct_data] = lcsv.load_fwd_csv_data(csv_file)
            subject = []

            # thr_data['thrtp_db'] = np.insert(thr_data['thrtp_db'], 0, np.NaN)  # put NaNs at ends of array
            # thr_data['thrtp_db'] = np.append(thr_data['thrtp_db'], np.NaN)
            thr_data['thrmp_db'] = np.array(thr_data['thrmp_db'])
            thr_data['thrtp_db'] = np.array(thr_data['thrtp_db'])

            # TODO this is for testing purposes only
            # thr_data['thrmp_db'] -= 5.0
            # thr_data['thrtp_db'] -= 5.0

        else:  # use threshold data from a subject
            rposvals = electrodes['rpos']
            thr_data['thrtp_db'] = np.insert(thr_data['thrtp_db'], 0, np.nan)  # put NaNs at ends of array
            thr_data['thrtp_db'] = np.append(thr_data['thrtp_db'], np.nan)

            # Calculate offset to get closer to what the model can produce
            mp_offset_db = np.nanmean(thr_data['thrmp_db']) - np.nanmean(mono_thr)
            tp_offset_db = np.nanmean(thr_data['thrtp_db']) - np.nanmean(tripol_thr)
            # Use the monopolar offset
            # overall_offset_db = mp_offset_db
            # overall_offset_db = 0
            # Use the mean of monopolar and tripolar offsets
            overall_offset_db = np.mean([mp_offset_db, tp_offset_db])
            # offset_mult = 1.0
            thr_data['thrmp_db'] -= overall_offset_db
            thr_data['thrtp_db'] -= overall_offset_db

        thresh_summary_full[scen, 0, 0, :] = thr_data['thrmp_db']
        thresh_summary_full[scen, 0, 1, :] = thr_data['thrtp_db']
        survvals = np.empty(NELEC)
        survvals[:] = np.nan
        ct_data = {'stdiameter': [], 'scala': [], 'elecdist': [], 'espace': espace, 'type': [], 'insrt_base': [],
                   'insert_apex': []}
        radius = 1.0
        ct_data['stdiameter'] = radius * 2.0 * (np.zeros(NELEC) + 1.0)

        rposvals = electrodes['rpos']  # save this for later  # TODO possibly not needed here
        saverposvals = rposvals

        cochlea_radius = ct_data['stdiameter'] / 2.0
        if np.isnan(thr_data['thrtp_sigma']) or thr_data['thrtp_sigma'] < 0.75 or thr_data['thrtp_sigma'] > 1:
            print('The sigma value for the TP configuration is invalid.')
            exit()

        fradvec = (ct_data['stdiameter'] / 2)  # smooth the radius data!!
        fradvec = sig.filtfilt(np.hanning(5) / sum(np.hanning(5)), 1, fradvec)
        simParams['cochlea']['radius'] = fradvec
        avec = np.arange(0, 1.005, 0.01)  # create the neuron count to neuron spikes transformation
        rlvec = NEURONS['coef'] * (np.power(avec, 2.0)) + (1 - NEURONS['coef']) * avec
        rlvec = NEURONS['neur_per_clust'] * np.power(rlvec, NEURONS['power'])
        rlvltable = np.stack((avec, rlvec))
        simParams['neurons']['rlvl'] = rlvltable

        # Construct the simParams structure
        simParams['cochlea'] = COCHLEA
        simParams['electrodes'] = electrodes
        simParams['channel'] = CHANNEL
        simParams['channel']['sigma'] = 0.0
        thresholds = np.empty(NELEC)  # Array to hold threshold data for different simulation values
        thresholds[:] = np.nan
        fitrposvals = np.zeros(NELEC)
        fitsurvvals = np.zeros(NELEC)
        par = Parameters()
        initvec = []

        if fit_mode == 'combined':  # Optimize survival and rpos to match MP and TP thresholds
            # Loop on electrodes, fitting rpos and survival fraction at each location
            nsols = np.zeros(NELEC)  # number of solutions found from the 2D maps

            for i in range(1, NELEC - 1):  # Fit params for each electrode (except ends, where there is no tripol value)
                mptarg = thr_data['thrmp_db'][i]
                tptarg = thr_data['thrtp_db'][i]

                fig4, ax5 = plt.subplots()
                ax5 = plt.contour(rpos_grid_vals, surv_grid_vals[2:], mono_thr[2:, :], [mptarg], colors='green')
                ax5.axes.set_xlabel('Rpos (mm)')
                ax5.axes.set_ylabel('Survival fraction')
                ax6 = plt.contour(rpos_grid_vals, surv_grid_vals[2:], tripol_thr[2:, :], [tptarg], colors='red')
                ax6.axes.set_xlabel('Rpos (mm)')
                ax6.axes.set_ylabel('Survival fraction')
                ax5.axes.xaxis.set_label('Rpos (mm)')
                ax5.axes.yaxis.set_label('Survival fraction')
                mpcontour = ax5.allsegs[0]
                tpcontour = ax6.allsegs[0]
                if not ifPlotGuessContours:
                    plt.close(fig4)

                nmp = len(mpcontour[0])
                ntp = len(tpcontour[0])
                mpx = np.zeros(nmp)
                mpy = np.zeros(nmp)
                tpx = np.zeros(ntp)
                tpy = np.zeros(ntp)

                for j in range(0, nmp):  # Should be able to do this without for loops
                    mpx[j] = mpcontour[0][j][0]
                    mpy[j] = mpcontour[0][j][1]

                for j in range(0, ntp):
                    tpx[j] = tpcontour[0][j][0]
                    tpy[j] = tpcontour[0][j][1]

                x, y = intsec.intersection(mpx, mpy, tpx, tpy)  # find intersection(s)

                # How many intersections? 0, 1 or more?  If single intersection use those values
                nsols[i] = len(x)
                if nsols[i] == 0:
                    # no solution. This shouldn't happen with known scenarios,
                    # since the forward model calculated threshold.
                    if use_fwd_model:
                        print('electrode: ', i, ';  no solution, but this is a known scenario')
                        print('min thr in 2D: ', np.min(mono_thr[2:, :]), ' and ', np.min(tripol_thr[2:, :]))

                        print('max thr in 2D: ', np.max(mono_thr[2:, :]), ' and ', np.max(tripol_thr[2:, :]))
                        print(' and goals are: ', mptarg, ' and ', tptarg)
                        # exit()

                    mp_idx, tp_idx = find_closest(mpx, mpy, tpx, tpy)

                    # rp_guess = mpx[mp_idx]
                    # sv_guess = mpy[mp_idx]
                    # rp_guess = tpx[tp_idx]  # Use tripolar best guess
                    # sv_guess = tpy[tp_idx]

                    rp_guess = (mpx[mp_idx] + tpx[tp_idx]) / 2.0
                    sv_guess = (mpy[mp_idx] + tpy[tp_idx]) / 2.0
                    ax_guess = plt.plot(rp_guess, sv_guess, 'x')

                    # # print("no solutions. Closest: ", mp_idx, ' and ', tp_idx,
                    #       ' , leading to guesses of (position, survival): ', rp_guess, sv_guess)

                    print('electrode: ', i, ";  no solutions. Closest: ", mp_idx, ' and ', tp_idx,
                          ' , OVERRIDING to: (position, survival): ', rp_guess, sv_guess)

                elif nsols[i] == 1:  # unique solution
                    print('electrode: ', i, ';  one solution: ', x, y)
                    rp_guess = x
                    sv_guess = y
                    ax_guess = plt.plot(rp_guess, sv_guess, 'x')

                else:  # multiple solutions
                    print('electrode: ', i, ';  ', nsols[i], ' solutions: ', x, ' and: ', y)
                    which_sols = np.zeros((4, int(nsols[i])))  # array for solutions and best fit
                    for sol in range(int(nsols[i])):  # Try all potential solutions; keep best
                        rp_guess = x[sol]
                        sv_guess = y[sol]

                        # Put variables into Parameters
                        # par.add('rpos_val', value=rp_guess, min=-0.8, max=0.8)
                        # par.add('surv_val', value=sv_guess, min=0.2, max=0.9)

                        # do fit, here with the default leastsq algorithm
                        # minner = Minimizer(objectivefunc_lmfit, par, fcn_args=(sigmaVals, simParams,
                        # fp, act_vals, thr_data, i))
                        # result = minner.minimize()
                        # which_sols[0, sol] = result.params['rpos_val']
                        # which_sols[1, sol] = result.params['surv_val']
                        # which_sols[2, sol] = np.mean(result.residual)
                        rposweight = 2
                        survweight = 0.0
                        if i > 1:  # calculate distance from previous coords
                            which_sols[3, sol] = np.sqrt(rposweight * (rp_guess - fitrposvals[i - 1]) ** 2 +
                                                         survweight * (sv_guess - fitsurvvals[i - 1]) ** 2)

                    # figure out which is best
                    print("which is best? use average")
                    # First attempt: pick rpos and surv that are closest to previous electrode (if there is one)
                    # min_val = np.amin(which_sols[3, :])
                    # closest_idx = np.where(which_sols[3, :] == min_val)[0]
                    # if len(closest_idx) > 1:
                    #     closest_idx = 0
                    # rp_guess = x[closest_idx]
                    # sv_guess = y[closest_idx]
                    # print('Closest is # ', closest_idx, ' , and guesses are: ', rp_guess, sv_guess)
                    # # print('And residual is: ', which_sols[2, closest_idx])
                    rp_guess = np.mean(x)
                    sv_guess = np.mean(y)

                    ax_guess = plt.plot(rp_guess, sv_guess, 'x')

                if ifPlotGuessContours:
                    plt.show()

                fitrposvals[i] = rp_guess
                fitsurvvals[i] = sv_guess

            # fix values for first and last electrodes as identical to their neighbors
            fitrposvals[0] = fitrposvals[1]
            fitrposvals[-1] = fitrposvals[-2]
            fitsurvvals[0] = fitsurvvals[1]
            fitsurvvals[-1] = fitsurvvals[-2]

            if not tp_extend:
                initvec = np.append(fitrposvals[1:-1], fitsurvvals[1:-1])
                for i, val in enumerate(initvec):  # place values in to the par object
                    if i < NELEC - 2:
                        #lb = -0.95  # lower and upper bounds for position
                        #ub = 0.95
                        lb = -0.95
                        ub = 0.95
                    elif NELEC - 2 <= i < 2 * (NELEC - 2):
                        lb = 0.0  # density
                        ub = 1.0
                    else:  # TODO leftover from tryng to fit external resistivity. Can remove
                        lb = 125
                        ub = 2500

                    par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)

                # end block for single electrode fits during the combined fit

            else:
                initvec = np.append(fitrposvals, fitsurvvals)

                for i, val in enumerate(initvec):  # place values in to the par object
                    if i < NELEC:
                        #  lb = -0.95  # lower and upper bounds for position
                        #   ub = 0.95
                        lb = -0.95
                        ub = 0.95
                    elif NELEC <= i < 2 * NELEC:
                        lb = 0.0  # density
                        ub = 1.0

                    par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)

            # end block for single electrode fits during the combined fit

        elif fit_mode == 'rpos':  #TODO may need fixing later on # fit rpos only; hold survival fixed as the values loaded from the scenario
            if use_minimizer == True:
                initvec = np.append(fitrposvals[1:-1], fitsurvvals[1:-1])  ## this line needs attention
                for i, val in enumerate(initvec):  # place values in to the par object
                    if i < NELEC - 2:
                        lb = -0.95  # lower and upper bounds for position
                        ub = 0.95
                    elif NELEC - 2 <= i < 2 * (NELEC - 2):
                        lb = 0.0  # density
                        ub = 1.0
                    else:  # TODO leftover from tryng to fit external resistivity. Can remove
                        lb = 125
                        ub = 2500

                    par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)

                # end block for single electrode fits during the combined fit

            else:
                initvec = np.append(np.zeros(NELEC), survvals)
                for i, val in enumerate(initvec):
                    if i < 16:
                        lb = -0.85  # lower and upper bounds
                        ub = 0.85
                        par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)
                    else:
                        par.add('v_%i' % i, value=initvec[i], vary=False)

        elif fit_mode == 'survival':  ## TODO handle tp_extend if needed
            if use_fwd_model:
                rposvals = electrodes['rpos']
                #survVals, electrodes['rpos']
            else:
                rposvals = subject_data.subj_ct_data(subject)
            if not tp_extend:
                initvec = np.append(rposvals[1:-1], (np.ones(NELEC - 2) * 0.5))
                if use_minimizer == True:
                    for i, val in enumerate(initvec):  # place values in to the par object
                        if i < NELEC - 2:
                            lb = initvec[
                                     i] - 0.2  #0.01# # lower and upper bounds for position #max of lb coordinates max(),-1<--radius of cylinder
                            ub = initvec[i] + 0.2
                            if lb <= -1 * radius + .05:
                                lb = -1 * radius + .05
                        #  if ub < 0:
                        #ub = .95

                        elif NELEC - 2 <= i < 2 * (NELEC - 2):
                            lb = 0.0  # density
                            ub = 1.0

                        par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)
                else:
                    # end block for single electrode fits during the combined fit
                    for i, val in enumerate(initvec):
                        if i < 16:
                            lb = rposvals[i] - 0.2  #0.01 # lower and upper bounds
                            ub = rposvals[i] + 0.2
                            if lb <= -1 * radius:
                                lb = -1 * radius + .05
                            if ub <= -1 * radius + .95:
                                ub = -1 * radius + .95
                        else:
                            lb = 0.1
                            ub = 1.0
                        par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)
            else:
                initvec = np.append(rposvals, (np.ones(NELEC) * 0.5))
                if use_minimizer == True:
                    for i, val in enumerate(initvec):  # place values in to the par object
                        if i < NELEC:
                            lb = initvec[i] - 0.2  #CHANGE BACK 0.01# lower and upper bounds for position
                            ub = initvec[i] + 0.2
                            if lb <= -1 * radius:
                                lb = -1 * radius + .05
                            if ub <= -1 * radius + .95:
                                ub = -1 * radius + .95
                        else:
                            lb = 0.0  # density
                            ub = 1.0

                        par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)
                else:
                    # end block for single electrode fits during the combined fit
                    for i, val in enumerate(initvec):
                        if i < 16:
                            lb = rposvals[i] - 0.2  #CHANGE BACK # lower and upper bounds
                            ub = rposvals[i] + 0.2
                            if lb <= -1 * radius:
                                lb = -1 * radius + .05
                            if ub <= -1 * radius + .95:
                                ub = -1 * radius + .95
                        else:
                            lb = 0.1
                            ub = 1.0
                        par.add('v_%i' % i, value=initvec[i], min=lb, max=ub)

        if use_minimizer:  # Now do the main fitting for all electrodes at once
            minner = Minimizer(objectivefunc_lmfit_all, par, nan_policy='omit',
                               fcn_args=(sigmaVals, simParams, fp, act_vals, thr_data))

            if use_fwd_model:
                result = minner.minimize(method='least_squares', ftol=fit_tol, diff_step=0.1)

            else:  # use CT data
                result = minner.minimize(method='least_squares', ftol=fit_tol, diff_step=0.1)

            if not tp_extend:  # store the results in the right place
                for i in range(NELEC - 2):
                    vname = 'v_%i' % i
                    fitrposvals[i + 1] = result.params[vname]
                    vname = 'v_%i' % (i + NELEC - 2)
                    fitsurvvals[i + 1] = result.params[vname]

                fitrposvals[0] = fitrposvals[1]  #potentially check this
                fitrposvals[-1] = fitrposvals[-2]
                fitsurvvals[0] = fitsurvvals[1]
                fitsurvvals[-1] = fitsurvvals[-2]


            else:
                for i in range(NELEC):
                    vname = 'v_%i' % i
                    fitrposvals[i] = result.params[vname]
                    vname = 'v_%i' % (i + NELEC)
                    fitsurvvals[i] = result.params[vname]

        else:  # use standard scipy.minimize
            bnds = ((-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9),
                    (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9),
                    (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (-0.9, 0.9), (0.05, 0.95), (0.05, 0.95),
                    (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95),
                    (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95), (0.05, 0.95),
                    (0.05, 0.95), (0.05, 0.95))

            result = opt.minimize(objectivefunc_minimize_all, initvec, args=(sigmaVals, simParams, fp, act_vals,
                                                                             thr_data), method='SLSQP', jac=None,
                                  bounds=bnds, options={'ftol': fit_tol})

            # store the results in the right place
            print('minimize finished. Message is: ', result.message)
            fitrposvals = result.x[0:NELEC]
            fitsurvvals = result.x[NELEC:]

        simParams['electrodes']['rpos'] = fitrposvals

        # report last fit
        if use_minimizer:
            report_fit(result)
            result.params.pretty_print()
        else:
            print('fitrposvals: ', fitrposvals)
            print('fitsurvvals: ', fitsurvvals)

        simParams['electrodes']['rpos'] = fitrposvals
        simParams['neurons']['nsurvival'] = surv_full.surv_full(simParams['electrodes']['zpos'], fitsurvvals,
                                                                simParams['grid']['z'])

        simParams['channel']['sigma'] = 0  # sigvals[0]
        simParams['channel']['behavior'] = 0
        thrsimmp = gt.get_thresholds(act_vals, fp, simParams)
        simParams['channel']['sigma'] = sigmaVals[1]
        #simParams['channel']['behavior'] = 1 for M levels
        simParams['channel']['behavior'] = 0
        thrsimtp = gt.get_thresholds(act_vals, fp, simParams)

        errvals = np.abs(np.subtract(thrsimmp[0], thr_data['thrmp_db'])), np.abs(np.subtract(thrsimtp[0][1:NELEC - 1],
                                                                                             thr_data['thrtp_db'][
                                                                                             1:NELEC - 1]))
        thrsim = [thrsimmp[0]], [thrsimtp[0]]
        thrtargs = [[thr_data['thrmp_db']], [thr_data['thrtp_db']]]

        # Summary data. Monopolar and tripolar errors saved separately.
        mean_mp_thr_err = np.nanmean(np.abs(np.array(thrsim[0]) - np.array(thrtargs[0])))
        mean_tp_thr_err = np.nanmean(np.abs(np.array(thrsim[1]) - np.array(thrtargs[1])))

        thresh_err_summary[scen, 0] = mean_mp_thr_err
        thresh_err_summary[scen, 1] = mean_tp_thr_err
        thresh_summary_full[scen, 1, 0, :] = np.array(thrsim[0])
        thresh_summary_full[scen, 1, 1, :] = np.array(thrsim[1])

        if use_fwd_model:
            [survvals, rposvals, espace] = s_scen.set_scenario(scenario, NELEC)
            if not tp_extend:
                rposerrs = np.subtract(rposvals[1:-1], fitrposvals[1:-1])
                survivalerrs = np.abs(np.subtract(survvals[1: -1], fitsurvvals[1: -1]))
            else:
                rposerrs = np.subtract(rposvals, fitrposvals)
                survivalerrs = np.abs(np.subtract(survvals, fitsurvvals))

            rpos_summary.append([rposvals, fitrposvals])
            rpos_err_metric = np.mean(np.abs(rposerrs))
            rpos_err_summary[scen] = rpos_err_metric
            surv_err_summary[scen] = np.mean(np.abs(survivalerrs))

            if not tp_extend:
                distvals = 1 - rposvals[1:-1]
                if np.std(distvals) == 0.0:
                    distvals += np.random.normal(0.0, 1e-20, size=len(distvals))  # to prevent a nan from stats.pearsonr
                    # distvals[-1] += 1e-8  # to prevent a nan from stats.pearsonr
                # [dist_corr[scen], dist_corr_p[scen]] = stats.pearsonr(1 - rposvals[1:-1], 1 - fitrposvals[1:-1])
                [dist_corr[scen], dist_corr_p[scen]] = stats.pearsonr(distvals, 1 - fitrposvals[1:-1])

            else:
                [dist_corr[scen], dist_corr_p[scen]] = stats.pearsonr(1 - rposvals, 1 - fitrposvals)

            # Save values in CSV format
            save_file_name = INVOUTPUTDIR + scenario + '_fitResults_' + fit_mode + '.csv'
        else:  # for subjects
            ct_vals = subject_data.subj_ct_data(subject)
            survivalerrs = np.empty(NELEC)

            if np.any(ct_vals):
                if not tp_extend:
                    rposerrs = np.abs(np.subtract(fitrposvals[1:-1], ct_vals[1:-1]))
                    rpos_summary.append([fitrposvals[1:-1], ct_vals[1:-1]])
                    rpos_err_metric = np.mean(np.abs(rposerrs))
                    rpos_err_summary[scen] = rpos_err_metric

                    [dist_corr[scen], dist_corr_p[scen]] = stats.pearsonr(1 - ct_vals[1:-1], 1 - fitrposvals[1:-1])

                else:
                    rposerrs = np.abs(np.subtract(fitrposvals, ct_vals))
                    rpos_err_metric = np.mean(np.abs(rposerrs))
                    rpos_summary.append([fitrposvals, ct_vals])
                    rpos_err_summary[scen] = rpos_err_metric
                    [dist_corr[scen], dist_corr_p[scen]] = stats.pearsonr(1 - ct_vals, 1 - fitrposvals)

            else:
                rposerrs = np.empty(NELEC)
                rpos_err_metric = np.nan
                rpos_err_summary[scen] = np.nan

            # Save values in CSV format
            save_file_name = INVOUTPUTDIR + subject + '_fitResults_' + fit_mode + '.csv'

        # Save the data for this scenario
        with open(save_file_name, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',', quotechar='"')
            if use_fwd_model:
                header = ['Electrode', 'ThreshMP', 'ThreshTP', 'Fitted ThreshMP', 'Fitted ThreshTP',
                          'Rposition', 'RpositionFit', 'Survival', 'Fitted Survival', 'RposError', 'SurvError',
                          'offset (NAN)']
            else:
                header = ['Electrode', 'ThreshMP', 'ThreshTP', 'Fitted ThreshMP', 'Fitted ThreshTP',
                          'CT_position', 'RpositionFit', 'Survival (Nan)', 'Fitted survival', 'RposError',
                          'SurvError (Nan)', 'offset']

            data_writer.writerow(header)
            for row in range(NELEC):
                t1 = row
                t2 = thr_data['thrmp_db'][row]
                t3 = thr_data['thrtp_db'][row]
                t4 = thrsimmp[0][row]
                t5 = thrsimtp[0][row]
                if use_fwd_model:
                    t6 = rposvals[row]
                else:
                    if np.any(ct_vals):
                        t6 = ct_vals[row]
                    else:
                        t6 = np.nan
                t7 = fitrposvals[row]
                if use_fwd_model:
                    t8 = survvals[row]
                else:
                    t8 = np.nan
                t9 = fitsurvvals[row]
                if row == 0 or row == NELEC - 1:
                    t10 = np.nan
                else:
                    t10 = rposerrs[row - 1]  # with tp_extend, only has 14 values
                if use_fwd_model:
                    if row == 0 or row == NELEC - 1:
                        t11 = np.nan
                        t12 = np.nan
                    else:
                        t11 = survivalerrs[row - 1]

                else:
                    t11 = np.nan
                    t12 = overall_offset_db

                data_writer.writerow([t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12])

            # Done with the values for each electrode
            thr_err_mp = np.abs(np.subtract(thrsim[0][0][:], thrtargs[0][0][:]))
            ovrl_mean_mp_err = np.mean(thr_err_mp)
            thr_err_tp = np.abs(np.subtract(thrsim[1][0][:], thrtargs[1][0][:]))
            ovrl_mean_tp_err = np.mean(thr_err_tp)

            data_writer.writerow(" ")  # blank line
            data_writer.writerow(['Mean: ', ' ', 'MP_err: ', ovrl_mean_mp_err, ' TP_err: ', ovrl_mean_tp_err])

        data_file.close()

        # Save values in npy format
        save_file_name = INVOUTPUTDIR + scenario + '_fitResults_' + fit_mode + '.npy'
        if use_fwd_model:
            np.save(save_file_name,
                    np.array([sigmaVals, rposvals, survvals, thrsim, thrtargs, initvec, [fitrposvals, fitsurvvals],
                              rposerrs, rpos_err_metric, survivalerrs], dtype=object))
        else:
            np.save(save_file_name,
                    np.array([sigmaVals, rposvals, survvals, thrsim, thrtargs, initvec, [fitrposvals, fitsurvvals],
                              rposerrs, rpos_err_metric, survivalerrs, ct_vals], dtype=object))

        np.save(INVOUTPUTDIR + scenario + '_simParams_inv', simParams)

        # Make plots
        if ifPlot:
            if use_fwd_model:
                txt_string = scenario
            else:
                txt_string = subject

            plot_inverse_results.plot_inverse_results(use_fwd_model, txt_string, unsupervised, fit_mode)
        if not if_save_npy:
            os.remove(save_file_name)
            os.remove(INVOUTPUTDIR + scenario + '_simParams_inv.npy')

    # Now we are done with the loop on scenarios/subjects
    pr.disable()  # stop the profiler
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(10)
    print(s.getvalue())

    print('Completed inverse model on ', num_scen, ' scenarios in total time (s) of: ', ps.total_tt)

    # save summary data in a CSV file
    # scenario, threshold error metric, and for CT cases: rpos error matrix
    # TODO: provide more info on inverse fit type in summary filename
    # Right now, each run of the inverse model overwrites the last summary. We should keep these results,
    # and somehow identify them as subject or scenario or something in the filename.
    print('saving summary data in CSV form')
    summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.csv'
    with open(summary_file_name, mode='w') as data_file:
        data_writer = csv.writer(data_file, delimiter=',', quotechar='"')
        header = ['Scenario/subject', 'Mean MP Threshold Error', 'Mean TP Threshold Error', 'Position Error',
                  'Density error', 'Dist correlation', 'Dist corr p']
        data_writer.writerow(header)
        avg_rpos_err = 0.0
        for row in range(num_scen):
            if is_scenario(scenarios[row]):
                data_writer.writerow([scenarios[row], '%.4f' % float(thresh_err_summary[row, 0]),
                                      '%.4f' % float(thresh_err_summary[row, 1]),
                                      '%.4f' % float(rpos_err_summary[row]),
                                      '%.4f' % float(surv_err_summary[row]),
                                      '%.4f' % float(dist_corr[row]),
                                      '%.5f' % float(dist_corr_p[row])])
            else:
                data_writer.writerow([scenarios[row], '%.4f' % float(thresh_err_summary[row, 0]),
                                      '%.4f' % float(thresh_err_summary[row, 1]),
                                      '%.4f' % float(rpos_err_summary[row]), '%.4f' % np.nan,
                                      '%.4f' % float(dist_corr[row]), '%.5f' % float(dist_corr_p[row])])

            avg_rpos_err += rpos_err_summary[row]

        avg_rpos_err /= num_scen  # convert from sum to average

        # Loop on scenarios/subjects done. Now provide summary data across all
        data_writer.writerow(' ')
        text = 'Summary. Mean rpos errorr: ' + '%.3f' % avg_rpos_err
        data_writer.writerow([text])
        data_file.close()

    # save summary binary data file
    # scenario, threshold error metric, and for CT cases: rpos error metric
    if if_save_npy:
        summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.npy'
        np.save(summary_file_name, np.array([scenarios, [thresh_summary_full], rpos_summary], dtype=object))

    print('Done with saving summary files')


if __name__ == '__main__':
    inverse_model_combined(mode='main')  # main or 'survey'
