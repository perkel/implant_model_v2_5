#  Fig_fit_summary.py
#  David Perkel 30 March 2024
import numpy as np

from common_params import *  # import common values across all models
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib as mpl
import seaborn as sns
import csv
import subject_data
import scipy.stats as stats
import pandas as pd

def is_scenario(scen):  # test whether this is a scenario or subject
    # if this scenario is a subject, set use_forward_model to be false
    if (scen[0] == 'A' or scen[0] == 'S') and scen[1:3].isnumeric():
        return False  # it's a subject
    else:
        return True

def read_inv_summary(res, na_scen):
    # Reads a summary file, and tests whether average rpos error is less than chance based on shuffling
    # You need to run the inverse model for all subjects before making this figure

    # construct correct path for this resistivity

    # INV_OUT_PRFIX = 'INV_OUTPUT/'
    # R_TEXT = 'R' + str(round(res))
    # INVOUTPUTDIR = INV_OUT_PRFIX + R_TEXT + ACTR_TEXT + STD_TEXT + TARG_TEXT
    # new_dir_suffix = 'R%d' % res + '_' + 'std_%.1f' % ACT_STDREL + '_thr_%d' % THRTARG + '/'
    new_dir_suffix = 'RE%d' %R_EXT + '_RI%d' % res + 'std_%.1f' % ACT_STDREL + '_thr_%d' % THRTARG + '/'

    INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix

    summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.npy'
    print(summary_file_name)
    [scenarios, thresh_summ_all, rpos_summary] = np.load(summary_file_name, allow_pickle=True)
    if na_scen > 0:
        scenarios = scenarios[na_scen:]  # Trim off the artificial scenarios, leaving just the subjects
    nscen = len(scenarios)
    n_elec = NELEC
    rpos_vals = np.zeros((nscen, n_elec))
    rpos_fit_vals = np.zeros((nscen, n_elec))
    thresh_err_summary = np.zeros((nscen, 2))
    rpos_err_summary = np.zeros(nscen)
    density_err_summary = np.zeros(nscen)
    dist_corr = np.zeros(nscen)
    dist_corr_p = np.zeros(nscen)

    for i, scen in enumerate(scenarios):
        print('scenario: ', scen)
        rpos_fit_vals[i, 1:-1] = rpos_summary[i+3][0]
        rpos_vals[i, 1:-1] = rpos_summary[i+3][1]

    # get detailed data from the CSV summary file
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
        return [np.asarray(thresh_summ_all[0]), thresh_err_summary, rpos_fit_vals[:, 1:-1], rpos_vals[:, 1:-1], rpos_err_summary,
                aaa_temp, dist_corr, dist_corr_p]
    else:
        return [np.asarray(thresh_summ_all[0]), thresh_err_summary, rpos_fit_vals, rpos_vals, rpos_err_summary, aaa_temp,
            dist_corr, dist_corr_p]

def fig9_summary():
    # Constants
    label_ypos = 1.05
    n_subj = 18
    nscen = len(scenarios)
    n_artscen = nscen-n_subj # number of artificial scenarios, which should not be plotted here
    n_elec = 16
    mpl.rcParams['font.family'] = 'Arial'

    # Color values (from Matlab plotting for Fig. 8)
    fig9_colors = np.zeros((n_subj, 3))
    fig9_colors[0, :] = [0, 0, 135]  # S22
    fig9_colors[1, :] = [0, 0, 193]  # S27
    fig9_colors[2, :] = [0, 0, 255]  # S29
    fig9_colors[3, :] = [0, 3, 255]  # S38
    fig9_colors[4, :] = [5, 73, 255]  # S40
    fig9_colors[5, :] = [14, 131, 254]  # S41
    fig9_colors[6, :] = [24, 192, 255]  # S42
    fig9_colors[7, :] = [34, 255, 255]  # S43
    fig9_colors[8, :] = [53, 255, 193]  # S46
    fig9_colors[9, :] = [91, 255, 136]  # S27
    fig9_colors[10, :] = [140, 255, 83]  # S49
    fig9_colors[11, :] = [194, 255, 39]  # S50
    fig9_colors[12, :] = [255, 255, 9]  # S52
    fig9_colors[13, :] = [254, 195, 10]  # S53
    fig9_colors[14, :] = [253, 135, 6]  # S54
    fig9_colors[15, :] = [252, 79, 5]  # S55
    fig9_colors[16, :] = [252, 24, 6]  # S56
    fig9_colors[17, :] = [252, 0, 0]  # S57
    fig9_colors /= 255.0

    # Need data from 2 resistivities
    r_vals = [70.0, 250.0]

    # Layout figure
    fig1, axs1 = plt.subplots(3, 2, figsize=(8, 8))
    fig1.tight_layout(pad=3)

    plt.figtext(0.02, 0.95, 'A', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.95, 'B', color='black', size=20, weight='bold')
    plt.figtext(0.02, 0.63, 'C', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.63, 'D', color='black', size=20, weight='bold')
    plt.figtext(0.02, 0.33, 'E', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.33, 'F', color='black', size=20, weight='bold')
    plt.figtext(0.255, 0.92, 'R' + str(round(r_vals[0])), color='black', size=16)
    plt.figtext(0.755, 0.92, 'R' + str(round(r_vals[1])), color='black', size=16)

    (thr_sum_all_0, thresh_err_summary_0, rpos_fit_vals_0, rpos_vals_0, rpos_err_summary_0, aaa_temp, dist_corr_0,
     dist_corr_p_0) = read_inv_summary(r_vals[0], n_artscen)
    (thr_sum_all_1, thresh_err_summary_1, rpos_fit_vals_1, rpos_vals_1, rpos_err_summary_1, aaa_temp, dist_corr_1,
     dist_corr_p_1) = read_inv_summary(r_vals[1], n_artscen)

    thresh_err_summary = np.zeros((2, len(scenarios), 2))
    thresh_err_summary[0, 3:, :] = thresh_err_summary_0
    thresh_err_summary[1, 3:, :] = thresh_err_summary_1
    #print('rposvals: ', rpos_vals[0], rpos_vals[1])
    rposerrs = np.zeros((2, nscen, nscen, n_elec))

    # color = iter(cm.rainbow(np.linspace(0, 1, n_subj)))
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel A
        x = np.asarray(thr_sum_all_0[idx + n_artscen, 0, 0, :])
        y = np.asarray(thr_sum_all_0[idx + n_artscen, 1, 0, :])
        axs1[0, 0].plot(x, y, '.', color=fig9_colors[idx, :])
        axs1[0, 0].set_xlabel('Measured monopolar threshold (dB)')
        axs1[0, 0].set_ylabel('Fit monopolar threshold (dB)')
        axs1[0, 0].spines['top'].set_visible(False)
        axs1[0, 0].spines['right'].set_visible(False)
        axs1[0, 0].set_xlim([30, 70])
        axs1[0, 0].set_ylim([30, 70])

    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel C
        x = np.asarray(thr_sum_all_0[idx + n_artscen, 0, 1, :])
        y = np.asarray(thr_sum_all_0[idx + n_artscen, 1, 1, :])
        axs1[1, 0].plot(x, y, '.', color=fig9_colors[idx, :])
        axs1[1, 0].set_xlabel('Measured tripolar threshold (dB)')
        axs1[1, 0].set_ylabel('Fit tripolar threshold (dB)')
        axs1[1, 0].spines['top'].set_visible(False)
        axs1[1, 0].spines['right'].set_visible(False)

    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel B
        x = np.asarray(thr_sum_all_1[idx + n_artscen, 0, 0, :])
        y = np.asarray(thr_sum_all_1[idx + n_artscen, 1, 0, :])
        axs1[0, 1].plot(x, y, '.', color=fig9_colors[idx, :])
        axs1[0, 1].set_xlabel('Measured monopolar threshold (dB)')
        axs1[0, 1].set_ylabel('Fit monopolar threshold (dB)')
        axs1[0, 1].spines['top'].set_visible(False)
        axs1[0, 1].spines['right'].set_visible(False)
        axs1[0, 1].set_xlim([30, 70])
        axs1[0, 1].set_ylim([30, 70])

    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel D
        x = np.asarray(thr_sum_all_1[idx + n_artscen, 0, 1, :])
        y = np.asarray(thr_sum_all_1[idx + n_artscen, 1, 1, :])
        axs1[1, 1].plot(x, y, '.', color=fig9_colors[idx, :])
        axs1[1, 1].set_xlabel('Measured tripolar threshold (dB)')
        axs1[1, 1].set_ylabel('Fit tripolar threshold (dB)')
        axs1[1, 1].spines['top'].set_visible(False)
        axs1[1, 1].spines['right'].set_visible(False)

    # Now distance data
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel E
        ct_dist = 1 - rpos_vals_0[idx]
        fit_dist = 1 - rpos_fit_vals_0[idx]
        retval = subject_data.subj_thr_data(scen)
        espace = retval[3]
        if espace == 0.85 or espace == 1.1:
            axs1[2, 0].plot(ct_dist, fit_dist, '.', color=fig9_colors[idx, :])
            axs1[2, 0].set_xlabel('Measured electrode distance (mm)')
            axs1[2, 0].set_ylabel('Fit electrode distance (mm)')
            axs1[2, 0].spines['top'].set_visible(False)
            axs1[2, 0].spines['right'].set_visible(False)
            axs1[2, 0].set_ylim(0, 2.0)
            [slope, intercept] = np.polyfit(ct_dist, fit_dist, 1)
            minx = np.min(ct_dist)
            maxx = np.max(ct_dist)
            axs1[2, 0].plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=fig9_colors[idx, :])  # plot line


    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel F
        ct_dist = 1 - rpos_vals_1[idx]
        fit_dist = 1 - rpos_fit_vals_1[idx]
        retval = subject_data.subj_thr_data(scen)
        espace = retval[3]
        if espace == 0.85 or espace == 1.1:
            axs1[2, 1].plot(ct_dist, fit_dist, '.', color=fig9_colors[idx, :])
            axs1[2, 1].set_xlabel('Measured electrode distance (mm)')
            axs1[2, 1].set_ylabel('Fit electrode distance (mm)')
            axs1[2, 1].spines['top'].set_visible(False)
            axs1[2, 1].spines['right'].set_visible(False)
            axs1[2, 1].set_ylim(0, 2.0)
            [slope, intercept] = np.polyfit(ct_dist, fit_dist, 1)
            minx = np.min(ct_dist)
            maxx = np.max(ct_dist)
            axs1[2, 1].plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=fig9_colors[idx, :])  # plot line


    # # Best fit line to the data
    coeffs = np.polyfit(1-rpos_vals_0.flatten(), 1-rpos_fit_vals_0.flatten(), 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*2.0)
    axs1[2, 0].plot([0, 2], [start_pt, end_pt], color='black')
    # # Now for the second resistivity
    coeffs = np.polyfit(1-rpos_vals_1.flatten(), 1-rpos_fit_vals_1.flatten(), 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*2.0)
    axs1[2, 1].plot([0, 2.0], [start_pt, end_pt], color='black')
     # axs1[2, 1].tick_params(
    #     axis='x',           # changes apply to the x-axis
    #     which='both',       # both major and minor ticks are affected
    #     bottom=False,       # ticks along the bottom edge are off
    #     top=False,          # ticks along the top edge are off
    #     labelbottom=False)  # labels along the bottom edge are off
    #
    # # statistics
    # # res = stats.linregress(x, y)

    # Save and display
    figname = 'Fig9_fit_summary.pdf'
    plt.savefig(figname, format='pdf', pad_inches=0.1)
    plt.show()

    plt.figure()
    dist_corr = np.zeros(nscen)
    dist_corr_p = np.zeros(nscen)

    signif = []
    near_signif = []
    not_signif = []

    for idx, corr in enumerate(dist_corr_0):
        if dist_corr_p_0[idx] > 0.1:
            not_signif.append(idx+3)
        elif dist_corr_p_0[idx] <= 0.1 and dist_corr_p_0[idx] > 0.05:
            near_signif.append(idx)
        else:
            signif.append(idx)


    sns.swarmplot(dist_corr_0)
    # sns.swarmplot(dist_corr_0[signif], color='black')
    # sns.swarmplot(dist_corr_0[near_signif], color='gray')
    # sns.swarmplot(dist_corr_0[not_signif], color='black', marker='$\circ$', s=6)
    plt.show()




if __name__ == '__main__':
    fig9_summary()
