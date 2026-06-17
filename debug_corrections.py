# script to check on algorithm for correcting electrode distance fits
import numpy as np
import matplotlib.pyplot as plt
from common_params import *  # import common values across all models
import matplotlib as mpl
import csv


def mean_error(arg, sl, intercept):
    return (arg * sl) + intercept


# takes as argument the fit distance and the average slope and intercept for errors
# from a population of subjects and returns a corrected fit distance, constrained to between 0.05 and 1.95 mm
def correct_error(f_d, sl, inter):
    err = (f_d * sl) + inter
    corrected_dist = f_d - err
    return corrected_dist, err


def read_inv_summary(res, na_scen):
    # Reads a summary file
    # You need to run the inverse model for all subjects before making this figure
    new_dir_suffix = 'RE%d' %R_EXT + '_RI%d' % res + 'std_%.1f' % ACT_STDREL + '_thr_%d' % THRTARG + '/'
    INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix

    summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.npy'
    print(summary_file_name)
    [scenarios, thresh_summ_all, rpos_summary] = np.load(summary_file_name, allow_pickle=True)
    if na_scen > 0:  # number of artificial scenarios
        scenarios = scenarios[na_scen:]  # Trim off the artificial scenarios, leaving just the subjects
        thresh_summ_all = thresh_summ_all[0][na_scen:][:][:][:]
        rpos_summary = rpos_summary[na_scen:]
    nscen = len(scenarios)
    n_elec = NELEC
    rpos_vals = np.zeros((nscen, n_elec-2))
    rpos_fit_vals = np.zeros((nscen, n_elec-2))
    thresh_err_summary = np.zeros((nscen, 2))
    rpos_err_summary = np.zeros(nscen)
    density_err_summary = np.zeros(nscen)
    dist_corr = np.zeros(nscen)
    dist_corr_p = np.zeros(nscen)

    for i, scen in enumerate(scenarios):  # loop on scenarios
        print('scenario: ', scen)
        # for single scenario
        rpos_fit_vals[i] = rpos_summary[i][0]
        # rpos_fit_vals[i, 1:-1] = rpos_summary[:]
        rpos_vals[i] = rpos_summary[i][1]

    # get detailed data from the CSV summary file
    aaa_temp = []
    summary_csv_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.csv'
    with open(summary_csv_file_name, mode='r') as data_file:
        entire_file = csv.reader(data_file, delimiter=',', quotechar='"')
        for row, row_data in enumerate(entire_file):
            if row < 4:  # skip header row and three scenarios
                pass
            else:
                if row < 22:
                    [_, thresh_err_summary[row-na_scen - 1, 0], thresh_err_summary[row-na_scen -1, 1],
                     rpos_err_summary[row-na_scen - 1], aaa_temp, dist_corr[row-na_scen - 1],
                     dist_corr_p[row-na_scen - 1]] = row_data
                else:
                    pass
                # Note aaa_temp is a placeholder for the density error, which is not used

        data_file.close()

    if not tp_extend:

        return [np.asarray(thresh_summ_all[0]), thresh_err_summary, rpos_fit_vals, rpos_vals[:],
                rpos_err_summary,
                aaa_temp, dist_corr, dist_corr_p]
        # return [np.asarray(thresh_summ_all[0]), thresh_err_summary, rpos_fit_vals[:, 1:-1], rpos_vals[:, 1:-1], rpos_err_summary,
        #         aaa_temp, dist_corr, dist_corr_p]
    else:
        return [np.asarray(thresh_summ_all[0]), thresh_err_summary, rpos_fit_vals, rpos_vals, rpos_err_summary, aaa_temp,
            dist_corr, dist_corr_p]



# make an array of 10 equally spaced numbers, corresponding to ct_dist
ct_dist = np.arange(10)*0.2  # zero through 9, just to get a spread

n_subj = 18
nscen = len(scenarios)
n_artscen = nscen-n_subj # number of artificial scenarios, which should not be plotted here
n_elec = 16
mpl.rcParams['font.family'] = 'Arial'
# Need data from 2 resistivities (leftover from Fig 9)
r_vals = [70.0, 250.0]

(thr_sum_all_0, thresh_err_summary_0, rpos_fit_vals_0, rpos_vals_0, rpos_err_summary_0, aaa_temp, dist_corr_0,
     dist_corr_p_0) = read_inv_summary(r_vals[0], n_artscen)  # read summary file with inverse model output
# assume some linear relationship
slope = 0.5
intercept = .3

fit_dist = (ct_dist * slope) + intercept  # set fake fit distances

# add noise
noise = np.random.rand(10)*0.2
fit_dist += noise

# overrride with real data from S subjects
# idx = 14
# ct_dist = 1 - rpos_vals_0[idx]
# fit_dist = 1 - rpos_fit_vals_0[idx]

# ct_dist = 1 - np.array((-0.19800,0.20000,0.32200,0.24200,-0.08900,-0.31500,-0.26200,-0.22400,-0.22000,-0.33600,-0.28300,-0.11800,0.08400,0.29700))
# fit_dist = 1- np.array((0.81510,0.87759,0.69908,0.77841,0.71341,0.85028,0.71947,0.73667,0.42570,0.21770,0.06412,-0.05898,-0.03440,0.28787))

# 2nd example from real data
# ct_dist = 1 - np.array((-0.64246,-0.50943,-0.67470,-0.65000,-0.64319,-0.59292,-0.65106,-0.57798,-0.40092,-0.28155,-0.22772,-0.14563,-0.20737,-0.31250))
# fit_dist = 1 - np.array((-0.03004,-0.08718,0.16489,0.41516,0.34516,0.62860,0.51414,0.10703,0.17881,0.45885,0.58992,0.65894,0.25761,0.07673))



err_dist = fit_dist - ct_dist  # calculate errors

fig, ax = plt.subplots(3, 2, figsize=(8, 10))
fig.tight_layout(pad=3)
ax[0, 0].plot(ct_dist, fit_dist, '.')  # Plot fit dist v. ct dist
ax[0, 0].set_xlabel('ct dist')
ax[0, 0].set_ylabel('fit_dist')
ax[0, 0].plot((0, 2), (0, 2), ':')
ax[0, 0].set_xlim(0, 2)
ax[0, 0].set_ylim(0, 2)

ax[2, 1].plot(fit_dist, ct_dist, '.')  # Plot fit dist v. ct dist
ax[2, 1].set_ylabel('ct dist')
ax[2, 1].set_xlabel('fit_dist')
ax[2, 1].plot((0, 2), (0, 2), ':')
ax[2, 1].set_xlim(0, 2)
ax[2, 1].set_ylim(0, 2)


ax[0, 1].plot(fit_dist, err_dist, '.')  # Plot error v. fit dist, as in paper
ax[0, 1].set_xlabel('ct dist')
ax[0, 1].set_ylabel('err_dist')

[err_sl, err_int] = np.polyfit(fit_dist, err_dist, 1)  # linear fit slope and intercept

# override slope and intercept using average from real data
# err_sl = 0.58
# err_int = -0.83
# err_sl = 0.347  # from Prism fit to all points together
# err_int = -0.66
ax[0, 1].plot((np.min(fit_dist), np.max(fit_dist)), ((np.min(fit_dist)* err_sl) + err_int, (np.max(fit_dist) * err_sl) + err_int))
ax[0, 1].set_xlabel('fit dist')
ax[0, 1].set_xlim(0, 2)
ax[0, 1].set_ylim(-1.5, 1.5)

dist_corr, err = correct_error(fit_dist, err_sl, err_int)  # Call function to calculate correction

err_2 = fit_dist*err_sl + err_int  # leftover code from before I created the function

ax[1, 0].plot(fit_dist, err, '.')
ax[1, 0].set_xlabel('fit dist')
ax[1, 0].set_ylabel('err')
ax[1, 0].set_xlim(0, 2)
ax[1, 0].set_ylim(-1.5, 1.5)

dist_corr_2 = fit_dist - err  # leftover code from before I created the function
ax[1, 1].plot(ct_dist, dist_corr, '.')
ax[1, 1].plot((0, 2), (0, 2), ':')
ax[1, 1].set_xlabel('ct dist')
ax[1, 1].set_ylabel('dist_corr')
ax[1, 1].set_xlim(0, 2)
ax[1, 1].set_ylim(0, 2)

err_after_corr = dist_corr - ct_dist
# error after correction
ax[2, 0].plot(ct_dist, err_after_corr, '.')
ax[2, 0].set_xlabel('ct dist')
ax[2, 0].set_ylabel('err_after_corr')

plt.show()
