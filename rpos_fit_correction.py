#  rpos_fit_correction.py
#  David Perkel 28 Feb 2026
# based on Fig9_v5.py from Perkel et al. paper

from common_params import *  # import common values across all models
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv

def is_scenario(scen):  # test whether this is a scenario or subject
    # if this scenario is a subject, set use_forward_model to be false
    if (scen[0] == 'A' or scen[0] == 'S') and scen[1:3].isnumeric():
        return False  # it's a subject
    else:
        return True

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


def fig9_summary():
    # Constants
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

    # Need data from 2 resistivities (leftover from Fig 9)
    r_vals = [70.0, 250.0]

    # Layout figure
    fig1, axs1 = plt.subplots(2, 2, figsize=(8, 8))
    fig1.tight_layout(pad=3)

    plt.figtext(0.02, 0.95, 'A', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.95, 'B', color='black', size=20, weight='bold')
    plt.figtext(0.02, 0.47, 'C', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.47, 'D', color='black', size=20, weight='bold')
    plt.figtext(0.255, 0.92, 'Original fits', color='black', size=16, horizontalalignment='center')
    plt.figtext(0.755, 0.92, 'Original fit errors', color='black', size=16, horizontalalignment='center')
    plt.figtext(0.255, 0.48, 'Corrected distances', color='black', size=16, horizontalalignment='center')
    plt.figtext(0.755, 0.48, 'Error after correction', color='black', size=16, horizontalalignment='center')

    (thr_sum_all_0, thresh_err_summary_0, rpos_fit_vals_0, rpos_vals_0, rpos_err_summary_0, aaa_temp, dist_corr_0,
     dist_corr_p_0) = read_inv_summary(r_vals[0], n_artscen)  # read summary file with inverse model output

    # thresh_err_summary = np.zeros((2, len(scenarios), 2))
    # thresh_err_summary[0, 3:, :] = thresh_err_summary_0
    err_dist = np.zeros((nscen, n_elec-2))

    # Now distance data
    case = 0  # for debugging
    summary_error = np.zeros((nscen-n_artscen, 9))  # for CT and fit dist, slope and intercept, correction
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel A
        ct_dist = 1 - rpos_vals_0[idx]
        fit_dist = 1 - rpos_fit_vals_0[idx]
        err_dist[idx, :] = fit_dist - ct_dist
        axs1[0, 0].plot(ct_dist, fit_dist, '.', color=fig9_colors[idx, :])

        [slope, intercept] = np.polyfit(ct_dist, fit_dist, 1)  # linear fit
        [err_slope, err_intercept] = np.polyfit(fit_dist, err_dist[idx], 1)
        # summary_error[idx, 0] = ct_dist
        # summary_error[idx, 1] = fit_dist
        summary_error[idx, 2] = slope
        summary_error[idx, 3] = intercept
        summary_error[idx, 7] = err_slope
        summary_error[idx, 8] = err_intercept
        minx = np.min(ct_dist)
        maxx = np.max(ct_dist)
        # if idx == case:
        #     axs1[0, 0].plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=fig9_colors[idx, :])  # plot line
        axs1[0, 1].plot(fit_dist, err_dist[idx], '.', color=fig9_colors[idx, :])  # plot error v. fit dist

    axs1[0, 0].set_xlabel('Measured electrode distance (mm)')  # Cosmetic labeling
    axs1[0, 0].set_ylabel('Fit electrode distance (mm)')
    axs1[0, 0].spines['top'].set_visible(False)
    axs1[0, 0].spines['right'].set_visible(False)
    axs1[0, 0].set_ylim(0, 2.0)
    axs1[0, 1].set_xlabel('Fit electrode distance (mm)')
    axs1[0, 1].set_ylabel('Fit error (mm)')
    axs1[0, 1].spines['top'].set_visible(False)
    axs1[0, 1].spines['right'].set_visible(False)
    axs1[0, 1].set_ylim(-1.2, 1.2)

    mean_slope = np.mean(summary_error[:, 2])
    mean_intercept = np.mean(summary_error[:, 3])
    mean_err_slope = np.mean(summary_error[:, 7])
    mean_err_intercept = np.mean(summary_error[:, 8])

    #  plot the  line with average slope and intercept
    axs1[0, 1].plot((0, 2), (mean_err_intercept, (2*mean_err_slope) + mean_err_intercept), '-')

    # Now plot corrected values
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel F
        min_dist = 0.05
        max_dist = 1.95
        ct_dist = 1 - rpos_vals_0[idx]
        fit_dist = 1 - rpos_fit_vals_0[idx]

        # call function to make the correction and return corrected distances and errors
        corr_fit_dist, mean_error = correct_error(fit_dist, mean_err_slope, mean_err_intercept)
        # mean_error = (fit_dist*mean_err_slope) + mean_err_intercept  # obsolete
        # corr_fit_dist = fit_dist - mean_error  # obsolete

        # impose limits if needed
        min_args = np.argwhere(corr_fit_dist < min_dist)
        corr_fit_dist[min_args] = min_dist
        max_args = np.argwhere(corr_fit_dist > max_dist)
        corr_fit_dist[max_args] = max_dist

        error_after_correc = corr_fit_dist - ct_dist

        axs1[1, 0].plot(ct_dist, corr_fit_dist, '.', color=fig9_colors[idx, :])

        [slope, intercept] = np.polyfit(ct_dist, corr_fit_dist, 1)  # fit corrected distance to a line
        summary_error[idx, 5] = slope
        summary_error[idx, 6] = intercept
        minx = np.min(ct_dist)
        maxx = np.max(ct_dist)
        axs1[1, 0].plot((minx, maxx), ((minx*slope) + intercept, (maxx*slope) + intercept), '-',
                        c=fig9_colors[idx, :])  # plot line

        axs1[1, 1].plot(fit_dist, error_after_correc, '.', color=fig9_colors[idx, :])
        if idx == case:
            axs1[1, 0].plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=fig9_colors[idx, :])  # plot line

    # clean up labeling for panel C
    axs1[1, 0].set_xlabel('Measured electrode distance (mm)')
    axs1[1, 0].set_ylabel('Corrected fit electrode distance (mm)')
    axs1[1, 0].spines['top'].set_visible(False)
    axs1[1, 0].spines['right'].set_visible(False)
    axs1[1, 0].set_xlim(0, 2.0)
    axs1[1, 0].set_ylim(0.0, 2.0)
    axs1[1, 1].set_xlabel('Fit electrode distance (mm)')
    axs1[1, 1].set_ylabel('Corrected fit electrode distance (mm)')

    np.savetxt(INVOUTPUTDIR + 'error_correct.csv', (summary_error[:-3]), delimiter=',')  # save output

    # # Best fit line to the data
    coeffs = np.polyfit(1-rpos_vals_0.flatten(), 1-rpos_fit_vals_0.flatten(), 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*2.0)
    axs1[0, 0].plot([0, 2], [start_pt, end_pt], color='black')
    axs1[1, 1].set_ylim(-1.2, 1.2)
    axs1[1, 1].spines['top'].set_visible(False)
    axs1[1, 1].spines['right'].set_visible(False)
    # # Now for the second resistivity
    coeffs = np.polyfit(ct_dist.flatten(), corr_fit_dist.flatten(), 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*2.0)
    axs1[1, 0].plot([0, 2.0], [start_pt, end_pt], color='black')

    # Save and display
    figname = 'Fig9_fit_summary.pdf'
    plt.savefig(figname, format='pdf', pad_inches=0.1)
    plt.show()

    fig2, axs2 = plt.subplots(3, 3, figsize=(8, 8))
    fig2.tight_layout(pad=3)  # just 3 test cases
    for idx in [6, 7, 8]:
        row = idx - 6
        ct_dist = 1 - rpos_vals_0[idx]
        fit_dist = 1 - rpos_fit_vals_0[idx]

        corr_fit_dist, mean_error = correct_error(fit_dist, mean_err_slope, mean_err_intercept)
        # mean_error = (fit_dist * mean_err_slope) + mean_err_intercept
        # corr_fit_dist = fit_dist - mean_error
        axs2[row, 0].plot(ct_dist, fit_dist, '.', color=fig9_colors[idx, :])
        axs2[row, 0].set_xlabel('ct dist')
        axs2[row, 0].set_ylabel('Fit electrode distance (mm)')
        axs2[row, 0].plot((0, 2), (0, 2), '--')
        axs2[row, 0].plot((0, 2), (mean_intercept, (2*mean_slope) + mean_intercept), ':')
        axs2[row, 0].set_xlim(0, 2.0)
        axs2[row, 0].set_ylim(0.0, 2.0)

        axs2[row, 1].plot(fit_dist, mean_error, '.', color=fig9_colors[idx, :])
        axs2[row, 1].plot((np.min(fit_dist), np.max(fit_dist)), (mean_err_intercept, (2*mean_err_slope) + mean_err_intercept), ':', color='black')
        axs2[row, 1].set_xlabel('fit dist')
        axs2[row, 1].set_ylabel('error (mm)')

        axs2[row, 2].plot(ct_dist, corr_fit_dist, '.', color=fig9_colors[idx, :])
        axs2[row, 2].set_xlabel('ct dist')
        axs2[row, 2].set_ylabel('Corrected electrode distance (mm)')
        axs2[row, 2].set_xlim(0, 2.0)
        axs2[row, 2].set_ylim(0.0, 2.0)

    plt.show()





if __name__ == '__main__':
    fig9_summary()
