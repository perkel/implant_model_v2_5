#  dist_fit_correction.py
#  David Perkel 15 April 2026
# based on Fig9_v5.py from Perkel et al. paper
# Now using an optimization approach and linear or binomial fits to minimize overall error

from common_params import *  # import common values across all models
import matplotlib.pyplot as plt
import matplotlib as mpl
from lmfit import Minimizer, Parameters, report_fit
import csv
import sys

def is_scenario(scen):  # test whether this is a scenario or subject
    # if this scenario is a subject, set use_forward_model to be false
    if (scen[0] == 'A' or scen[0] == 'S') and scen[1:3].isnumeric():
        return False  # it's a subject
    else:
        return True

# takes as argument the fit distance and the average slope and intercept for errors
# from a population of subjects and returns a corrected fit distance, constrained to be between 0.05 and 1.95
def correct_error(f_d, vals):
    n_pts = len(vals)
    min_dist = 0.05
    max_dist = 1.95

    v = np.zeros(n_pts)
    for idx, key in enumerate(vals.keys()):
        v[idx] = vals[key]

    if n_pts == 2:
        corrected_dist = (f_d - v[1])/v[0]
        c_d2 = np.polyval(np.roots(v), f_d)
    elif n_pts == 3:
        corrected_dist = np.polyval(np.roots(v), f_d)
    else:
        # raise error
        raise RuntimeError('unknown # of parameters in correction')

    corrected_dist[corrected_dist < min_dist] = min_dist
    corrected_dist[corrected_dist > max_dist] = max_dist
    return corrected_dist

# takes as argument the fit distance and the average slope and intercept for errors
# from a population of subjects and returns a corrected fit distance, constrained to between 0.05 and 1.95 mm
def correct_error_lin(f_d, sl, inter):
    aaa = (f_d - inter)/sl

    return aaa


def read_inv_summary(res, na_scen):
    # Reads a summary file
    # You need to run the inverse model for all subjects before making this figure
    new_dir_suffix = 'RE%d' %R_EXT + '_RI%d' % res + 'std_%.1f' % ACT_STDREL + '_thr_%d' % THRTARG + '/'
    INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix

    summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results_Seattle.npy'
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
        if len(rpos_summary[i][0]) > 14:
            rpos_fit_vals[i] = rpos_summary[i][0][1:-1]
            rpos_vals[i] = rpos_summary[i][1][1:-1]
        else:
            rpos_fit_vals[i] = rpos_summary[i][0]

            # rpos_fit_vals[i, 1:-1] = rpos_summary[:]
            rpos_vals[i] = rpos_summary[i][1]

    # get detailed data from the CSV summary file
    aaa_temp = []
    summary_csv_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results_Seattle.csv'
    with open(summary_csv_file_name, mode='r') as data_file:
        entire_file = csv.reader(data_file, delimiter=',', quotechar='"')
        for row, row_data in enumerate(entire_file):
            if row < 4:  # skip header row and three scenarios
                pass
            else:
                if row < 20:
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


def cost_fcn(par, ct, fit):
    vals = par.valuesdict()
    nsubj = ct.shape[0]
    # extract params
    error = 0.0

    for i in range(nsubj):# loop on subjects
        # make correction
        corr_dist = correct_error(fit[i, :], vals)

        # compare ct to corrected fit
        error += np.mean(np.abs(np.subtract(corr_dist, ct[i, :])))

    error /= nsubj
    print('slope: ', vals['v_0'], ' and int: ', vals['v_1'], ' returning error: ', error)
    return error


import csv


def extract_nth_column(filepath, col: int = 1, skip_rows=1):
    """
    Extract the cot column from a CSV file.

    Args:
        filepath:  Path to the CSV file
        skip_rows: Number of header rows to skip (default: 1)

    Returns:
        List of values from the 4th column
    """
    results = []

    with open(filepath, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)

        # Skip the first n rows
        for _ in range(skip_rows):
            next(reader, None)

        for row in reader:
            if len(row) >= col:  # guard against short/empty rows
                results.append(float(row[col-1]))  # 0-indexed, so index 3 = 4th column

    return results


def extract_rpos(path, col, kip_rows=1):

    filepath = sys.argv[1] if len(sys.argv) > 1 else "data.csv"
    skip_rows = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    column = np.array(extract_nth_column(path, col, skip_rows))

    print(f"Extracted {len(column)} values from column " + str(col) + ":")
    for value in column:
        print(value)

    return column

def fig9_summary():
    # Constants
    n_subj = 18
    nscen = len(scenarios)
    n_artscen = nscen-n_subj # number of artificial scenarios, which should not be plotted here
    n_elec = 16
    mode = 'linear'
    # mode = 'quadratic' # alternative is linear
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
    plt.figtext(0.255, 0.92, 'Original fits', color='black', size=14, horizontalalignment='center')
    # plt.figtext(0.755, 0.92, 'Original fit errors', color='black', size=16, horizontalalignment='center')
    plt.figtext(0.255, 0.48, 'Optimized corrected distances', color='black', size=14, horizontalalignment='center')
    # plt.figtext(0.755, 0.48, 'Error after correction', color='black', size=16, horizontalalignment='center')

    (thr_sum_all_0, thresh_err_summary_0, rpos_fit_vals_0, rpos_vals_0, rpos_err_summary_0, aaa_temp, dist_corr_0,
     dist_corr_p_0) = read_inv_summary(r_vals[0], n_artscen)  # read summary file with inverse model output

    ct_dist = np.zeros((nscen-n_artscen, n_elec-2))
    fit_dist = np.zeros((nscen-n_artscen, n_elec-2))
    lincoefs = np.zeros((nscen-n_artscen, 2))
    binomcoeffs = np.zeros((nscen-n_artscen, 3))

    # Convert to distance and consolidate for all scenarios
    summary_error = np.zeros((nscen-n_artscen, 9))  # for CT and fit dist, slope and intercept, correction
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel A
        ct_dist[idx, :] = 1 - rpos_vals_0[idx]
        fit_dist[idx, :] = 1 - rpos_fit_vals_0[idx]
        axs1[0, 0].plot(ct_dist[idx, :], fit_dist[idx, :], '.', color=fig9_colors[idx, :], linewidth=2)
        axs1[0, 0].plot(ct_dist[idx, :], fit_dist[idx, :], '.', color=fig9_colors[idx, :], linewidth=2)

        lincoefs[idx, :] = np.polyfit(ct_dist[idx], fit_dist[idx], 1)  # linear fit
        binomcoeffs[idx, :] = np.polyfit(ct_dist[idx], fit_dist[idx], 2)

    # calculate fit of error v. fit dist from all points together
    sl_all, int_all = np.polyfit(ct_dist.flatten(), fit_dist.flatten(), 1)
    binomcoeffs_all = np.polyfit(ct_dist.flatten(), fit_dist.flatten(), 2)
    par = Parameters()
    initvec = []

    if mode == 'linear':
        par.add('v_0', value=sl_all, min=0.05, max=1.95)
        par.add('v_1', value=int_all, min=0.05, max=1.95)
    else:  # quadratic
        par.add('v_0', value=binomcoeffs_all[0], min=0.05, max=1.95)
        par.add('v_1', value=binomcoeffs_all[1], min=0.05, max=1.95)
        par.add('v_2', value=binomcoeffs_all[2], min=0.05, max=1.95)

    minner = Minimizer(cost_fcn, par, nan_policy='omit', fcn_args=(ct_dist, fit_dist))
    result = minner.minimize(method='least_squares', diff_step=0.1)

    result_sing = result

    print('optimal slope: ', result.params['v_0'].value, ' and intercept: ', result.params['v_1'].value)
    axs1[0, 0].set_xlabel('CT electrode distance (mm)')  # Cosmetic labeling
    axs1[0, 0].set_ylabel('Fit electrode distance (mm)')
    axs1[0, 0].spines['top'].set_visible(False)
    axs1[0, 0].spines['right'].set_visible(False)
    axs1[0, 0].set_ylim(0, 2.0)
    # unity line
    axs1[0, 0].plot((0, 2), (0, 2), ':')

    # Linear and quadratic fits
    minct = np.min(ct_dist.flatten())
    maxct = np.max(ct_dist.flatten())
    axs1[0, 0].plot((minct, maxct), np.polyval((sl_all, int_all), (minct, maxct)), '-', color='black')
    ct_range = np.linspace(minct, maxct, 50)
    pv = np.polyval(binomcoeffs_all, ct_range)
    axs1[0, 0].plot(ct_range, pv, '-')


    # axs1[0, 1].set_xlabel('Fit electrode distance (mm)')
    # axs1[0, 1].set_ylabel('Fit error (mm)')
    # axs1[0, 1].spines['top'].set_visible(False)
    # axs1[0, 1].spines['right'].set_visible(False)
    # axs1[0, 1].set_ylim(-2, 2)

    mean_slope = np.mean(summary_error[:, 2])
    mean_intercept = np.mean(summary_error[:, 3])
    mean_err_slope = np.mean(summary_error[:, 7])
    mean_err_intercept = np.mean(summary_error[:, 8])

    #  plot the  line with average slope and intercept
    # axs1[0, 1].plot((0, 2), (mean_err_intercept, (2*mean_err_slope) + mean_err_intercept), '-')

    corr_fit_dist = np.zeros((nscen-n_artscen, n_elec-2))
    corr_fit_dist_all = np.zeros((nscen-n_artscen, n_elec-2))
    # Now plot corrected values
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel F
        # call function to make the correction and return corrected distances and errors
        corr_fit_dist[idx, :] = correct_error(fit_dist[idx, :], result.params)
        axs1[1, 0].plot(ct_dist[idx, :], corr_fit_dist[idx, :], '.', color=fig9_colors[idx, :])

        vals_all = []
        if mode == 'linear':
            vals_all = {'v_0': sl_all, 'v_1': int_all}
        corr_fit_dist_all[idx, :] = correct_error(fit_dist[idx, :], vals_all)
        axs1[1, 1].plot(ct_dist[idx, :], corr_fit_dist_all[idx, :], '.', color=fig9_colors[idx, :])

        [slope, intercept] = np.polyfit(ct_dist[idx, :], corr_fit_dist[idx, :], 1)  # fit corrected distance to a line
        summary_error[idx, 5] = slope
        summary_error[idx, 6] = intercept
        minx = np.min(ct_dist[idx, :])
        maxx = np.max(ct_dist[idx, :])
        axs1[1, 0].plot((minx, maxx), ((minx*slope) + intercept, (maxx*slope) + intercept), '-',
                        c=fig9_colors[idx, :])  # plot line

        [slope, intercept] = np.polyfit(ct_dist[idx, :], corr_fit_dist_all[idx, :], 1)  # fit corrected distance to a line
        axs1[1, 1].plot((minx, maxx), ((minx * slope) + intercept, (maxx * slope) + intercept), '-',
                        c=fig9_colors[idx, :])  # plot line
    # clean up labeling for panel C
    axs1[1, 0].set_xlabel('Measured electrode distance (mm)')
    axs1[1, 0].set_ylabel('Corrected fit electrode distance (mm)')
    axs1[1, 0].spines['top'].set_visible(False)
    axs1[1, 0].spines['right'].set_visible(False)
    axs1[1, 0].set_xlim(0, 2.0)
    axs1[1, 0].set_ylim(0.0, 2.0)
    axs1[1, 0].plot((0, 2), (0, 2), ':')

    # clean up labeling for panel D
    axs1[1, 1].set_xlabel('Measured electrode distance (mm)')
    axs1[1, 1].set_ylabel('Corrected fit electrode distance (mm)')
    axs1[1, 1].set_title('Correction from all pts')
    axs1[1, 1].spines['top'].set_visible(False)
    axs1[1, 1].spines['right'].set_visible(False)
    axs1[1, 1].set_xlim(0, 2.0)
    axs1[1, 1].set_ylim(0.0, 2.0)
    axs1[1, 1].plot((0, 2), (0, 2), ':')

    np.savetxt(INVOUTPUTDIR + 'error_correct.csv', (summary_error[:-3]), delimiter=',')  # save output

    # # Best fit line to the data
    coeffs = np.polyfit(1-rpos_vals_0.flatten(), 1-rpos_fit_vals_0.flatten(), 1)
    # start_pt = coeffs[1]
    # end_pt = coeffs[1] + (coeffs[0]*2.0)
    # axs1[0, 0].plot([0, 2], [start_pt, end_pt], color='black')
    # # # Now for the second resistivity
    # coeffs = np.polyfit(ct_dist.flatten(), corr_fit_dist.flatten(), 1)
    # start_pt = coeffs[1]
    # end_pt = coeffs[1] + (coeffs[0]*2.0)
    # axs1[1, 0].plot([0, 2.0], [start_pt, end_pt], color='black')

    # Save and display
    figname = 'Fig9_fit_summary.pdf'
    plt.savefig(figname, format='pdf', pad_inches=0.1)
    plt.show()

    # Layout figure
    nr = 6
    nc = 3
    fig2, axs2 = plt.subplots(nr, nc, figsize=(7, 12))
    fig2.tight_layout(pad=3)
    for idx, scen in enumerate(scenarios[n_artscen:]):  # Panel F
        # print(idx)
        # call function to make the correction and return corrected distances and errors

        # corr_fit_dist[idx, :] = correct_error(fit_dist[idx, :], mean_slope, mean_intercept)
        # corr_fit_dist_all[idx, :] = correct_error(fit_dist[idx, :], result.params['v_0'].value, result.params['v_0'].value)
        row, col = divmod(idx, nc)
        # print(row, col)
        # axs2[row, col].plot(ct_dist, fit_dist, '.', color=fig9_colors[idx, :])
        axs2[row, col].plot(ct_dist[idx, :], fit_dist[idx, :], '.', color='red')
        coeffs = np.polyfit(ct_dist[idx, :], fit_dist[idx, :], 2)
        minct = np.min(ct_dist[idx, :])
        maxct = np.max(ct_dist[idx, :])
        ct_range = np.linspace(minct, maxct, 50)
        pv = np.polyval(coeffs, ct_range)
        # axs2[row, col].plot(ct_dist[idx, :], corr_fit_dist[idx, :], '.', color='blue')
        axs2[row, col].plot(ct_range, pv, '-', color=fig9_colors[idx, :])
        if idx == 0:
            axs2[row, col].set_title(scen + 'red=orig; blue=corrected')
        else:
            axs2[row, col].set_title(scen)
        axs2[row, col].set_xlabel('CT distance (mm)')
        axs2[row, col].set_ylabel('Fit distance (mm)')
        axs2[row, col].spines['top'].set_visible(False)
        axs2[row, col].spines['right'].set_visible(False)
        axs2[row, col].set_xlim(0, 2.0)
        axs2[row, col].set_ylim(0.0, 2.0)
        axs2[row, col].plot((0, 2), (0, 2), ':')

    plt.show()

    scen_single = "A022R"

    ## TODO
    col_to_extract = 7
    fp = INVOUTPUTDIR + scen_single + "_fitResults_combined.csv"
    # get correct column from INV_OUT CSV file and put the values in rpos_single
    rpos_single = extract_rpos(fp, col_to_extract, 1)
    # third figure for single case
    fig3 = plt.subplots(2, 2, figsize=(7, 12))
    # fig3.tight_layout(pad=3)


    dist_single = 1- rpos_single
    corr_fit_dist_single = correct_error(dist_single,result_sing.params)
    print(dist_single)
    print(corr_fit_dist_single)
    corr_rpos = 1 - corr_fit_dist_single

    # TODO this printing should make it easy to copy and past into subject_data
    for pos in corr_rpos:
        print(str(pos) + ', ')

    # Save to a corrected CSV file with "_corrected" appended to the file name





if __name__ == '__main__':
    fig9_summary()
