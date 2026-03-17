#  Fig9_v4.py
#  David Perkel 30 March 2024
# Updates 7 May 2024
import numpy as np

from common_params import *  # import common values across all models
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib as mpl
import seaborn as sns
import csv
import scipy.stats as stats
import pandas as pd


def read_inv_summary(res):
    # Reads a summary file, and tests whether average rpos error is less than chance based on shuffling
    # You need to run the inverse model for all subjects before making this figure

    # construct correct path for this resistivity
    R_TEXT = 'R' + str(round(res))
    INV_OUT_PRFIX = 'INV_OUTPUT/'
    INVOUTPUTDIR = INV_OUT_PRFIX + R_TEXT + ACTR_TEXT + STD_TEXT + TARG_TEXT

    summary_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.npy'

    [scenarios, rpos_summary] = np.load(summary_file_name, allow_pickle=True)
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
        rpos_fit_vals[i, :] = rpos_summary[i][0]
        rpos_vals[i, :] = rpos_summary[i][1]

    # get detailed data from the CSV summary file
    summary_csv_file_name = INVOUTPUTDIR + 'summary_inverse_fit_results.csv'
    with open(summary_csv_file_name, mode='r') as data_file:
        entire_file = csv.reader(data_file, delimiter=',', quotechar='"')
        for row, row_data in enumerate(entire_file):
            if row == 0:  # skip header row
                pass
            else:
                [_, thresh_err_summary[row-1, 0], thresh_err_summary[row-1, 1],
                 rpos_err_summary[row-1], aaa_temp, dist_corr[row-1], dist_corr_p[row-1]] = row_data
                # Note aaa_temp is a placeholder for the density error, which is not used

        data_file.close()

    return [thresh_err_summary, rpos_fit_vals, rpos_vals, rpos_err_summary, aaa_temp, dist_corr, dist_corr_p]

def fig9_summary():
    # Constants
    label_ypos = 1.05
    n_subj = 18
    nscen = len(scenarios)
    n_elec = 16
    mpl.rcParams['font.family'] = 'Arial'

    # Need data from 2 resistivities
    r_vals = [70.0, 375.0]

    # Layout figure
    
    fig1= plt.figure(figsize=(8, 5))
    gs = fig1.add_gridspec(20, 9)
    # fig1.tight_layout(pad=3)

    plt.figtext(0.06, 0.90, 'A', color='black', size=20, weight='bold')
    plt.figtext(0.65, 0.90, 'B', color='black', size=20, weight='bold')
    plt.figtext(0.06, 0.45, 'C', color='black', size=20, weight='bold')
    plt.figtext(0.5, 0.45, 'D', color='black', size=20, weight='bold')


    # get data
    thresh_err_summary_0, rpos_fit_vals_0, rpos_vals_0, rpos_err_summary_0, aaa_temp, dist_corr_0, dist_corr_p_0 =(
        read_inv_summary(r_vals[0]))
    thresh_err_summary_1, rpos_fit_vals_1, rpos_vals_1, rpos_err_summary_1, aaa_temp, dist_corr_1, dist_corr_p_1 =(
        read_inv_summary(r_vals[1]))

    thresh_err_summary = np.zeros((2, len(scenarios), 2))
    thresh_err_summary[0, :, :] = thresh_err_summary_0
    thresh_err_summary[1, :, :] = thresh_err_summary_1
    rpos_fit_vals = np.zeros((2, len(scenarios), n_elec))
    rpos_fit_vals[0, :] = rpos_fit_vals_0
    rpos_fit_vals[1, :] = rpos_fit_vals_1
    rpos_vals = np.zeros((2, nscen, n_elec))
    rpos_vals[0, :] = rpos_vals_0
    rpos_vals[1, :] = rpos_vals_1
    #print('rposvals: ', rpos_vals[0], rpos_vals[1])
    rposerrs = np.zeros((2, nscen, nscen, n_elec))
    dist_corr = np.zeros((2, nscen))
    dist_corr[0, :] = dist_corr_0
    dist_corr[1, :] = dist_corr_1
    dist_corr_p = np.zeros((2, nscen))
    dist_corr_p[0, :] = dist_corr_p_0
    dist_corr_p[1, :] = dist_corr_p_1

    # Loop on scenarios again. Compute pairwise mean absolute position error
    mean_errs = np.zeros([2, nscen, nscen])
    corr_vals = np.zeros([2, nscen, nscen])
    corr_p = np.zeros([2, nscen, nscen])

    axs0 = fig1.add_subplot(gs[0:9, 0:6])
    axs1 = fig1.add_subplot(gs[0:9, 7:])
    axs2 = fig1.add_subplot(gs[11:, 0:4])
    axs3 = fig1.add_subplot(gs[11:, 5:])
    axs0.plot(thresh_err_summary[0, :, 0], thresh_err_summary[0, :, 1], 'o', color='black', label="R 70")
    axs0.plot(thresh_err_summary[1, :, 0], thresh_err_summary[1, :, 1], 'o', color='blue', label="R 375")
    axs0.plot([0, 0.8], [0, 0.8], linestyle='dashed', color='black', label='y = x')  # line of slope 1.0
    axs0.set_xlabel('Monopolar threshold error (dB)')
    axs0.set_ylabel('Tripolar thresh. error (dB)')
    axs0.spines['top'].set_visible(False)
    axs0.spines['right'].set_visible(False)
    axs0.set_xlim([0, 1.5])
    axs0.set_ylim([-0.05, 0.85])
    axs0.legend()

    # Best fit line to the data
    coeffs = np.polyfit(thresh_err_summary[0, :, 0], thresh_err_summary[0, :, 1], 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*1.2)
    axs0.plot([0, 1.2], [start_pt, end_pt], color='black')
    # Now for the second resistivity
    coeffs = np.polyfit(thresh_err_summary[1, :, 0], thresh_err_summary[1, :, 1], 1)
    start_pt = coeffs[1]
    end_pt = coeffs[1] + (coeffs[0]*1.5)
    axs0.plot([0, 1.5], [start_pt, end_pt], color='blue')

    pos_err = pd.DataFrame(
        {
            "R70": rpos_err_summary_0[3:],
            "R375": rpos_err_summary_1[3:]
        }
    )

    sns.swarmplot(data=pos_err, ax=axs1, palette=['black', 'blue'])
    sns.lineplot(data=pos_err, ax=axs1)
    axs1.set_ylabel('Electrode distance error (mm)')
    axs1.spines['top'].set_visible(False)
    axs1.spines['right'].set_visible(False)

    # # plot median and interquartile range
    quartile = np.quantile(rpos_err_summary_0, [0.25, 0.5, 0.75])
    axs1.plot([-0.3, 0.3], [quartile[1], quartile[1]], 'k')  # horizontal line to indicate median
    axs1.plot([-0.15, 0.15], [quartile[0], quartile[0]], 'k')
    axs1.plot([-0.15, 0.15], [quartile[2], quartile[2]], 'k')
    axs1.plot([0, 0], [quartile[0], quartile[2]], 'k')  # vertical line
    axs1.set_xlim([-0.5, 1.5])

    quartile = np.quantile(rpos_err_summary_1, [0.25, 0.5, 0.75])
    axs1.plot([0.7, 1.3], [quartile[1], quartile[1]], 'b')  # horizontal line to indicate median
    axs1.plot([0.85, 1.15], [quartile[0], quartile[0]], 'b')
    axs1.plot([0.85, 1.15], [quartile[2], quartile[2]], 'b')
    axs1.plot([1, 1], [quartile[0], quartile[2]], 'b')  # vertical line
    axs1.set_xlim([-0.5, 1.5])

    axs1.set_xticks([0, 1], ['R 70', 'R 375'])

    axs2.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=True,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        labelbottom=True)  # labels along the bottom edge are off

    color = iter(cm.rainbow(np.linspace(0, 1, n_subj)))

    for subj in range(3, n_subj):
            c = next(color)
            ct_dist = 1 - rpos_vals[0, subj]
            fit_dist = 1 - rpos_fit_vals[0, subj]
            axs2.plot(ct_dist, fit_dist, '.', c=c)  # plot data points
            [slope, intercept] = np.polyfit(ct_dist, fit_dist, 1)
            minx = np.min(ct_dist)
            maxx = np.max(ct_dist)
            axs2.plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=c)  # plot line

    axs2.set_xlabel('CT estimate of electrode distance (mm)')
    axs2.set_ylabel('Fit electrode distance (mm)')
    axs2.set_ylim([0.0, 1.9])
    axs2.spines['top'].set_visible(False)
    axs2.spines['right'].set_visible(False)
    axs2.text(0.15, 1.7, "R 70")

    # best fit line for all points together
    x = 1 - rpos_vals[0, :, :].flatten()
    y = 1 - rpos_fit_vals[0, :, :].flatten()
    [slope, intercept] = np.polyfit(x, y, 1)
    start_pt, end_pt = axs2.get_xlim()
    axs2.plot([start_pt, end_pt], [(start_pt*slope) + intercept, (end_pt*slope) + intercept], 'black')
    print('For overall fit for R70, slope = %.2f' % slope, ' and intercept = %.2f' % intercept)

    # Now a similar scatter plot for the R70 data
    # Get data
    # Construct plot
    axs2.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=True,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        labelbottom=True)  # labels along the bottom edge are off

    color = iter(cm.rainbow(np.linspace(0, 1, n_subj)))

    for subj in range(3, n_subj):
        c = next(color)
        ct_dist = 1 - rpos_vals[1, subj]
        fit_dist = 1 - rpos_fit_vals[1, subj]
        axs3.plot(ct_dist, fit_dist, '.', c=c)  # plot data points
        [slope, intercept] = np.polyfit(ct_dist, fit_dist, 1)
        minx = np.min(ct_dist)
        maxx = np.max(ct_dist)
        axs3.plot((minx, maxx), (minx*slope + intercept, maxx*slope + intercept), '-', c=c)  # plot line

    axs3.set_xlabel('CT estimate of electrode distance (mm)')
    axs3.set_ylabel('Fit electrode distance (mm)')
    axs3.spines['top'].set_visible(False)
    axs3.spines['right'].set_visible(False)
    axs3.text(0.15, 1.7, "R 375")
    axs3.set_ylim([0.0, 1.9])

    # best fit line for all points together
    x = 1 - rpos_vals[1, 3:].flatten()
    y = 1 - rpos_fit_vals[1, 3:].flatten()
    [slope, intercept] = np.polyfit(x, y, 1)
    start_pt, end_pt = axs3.get_xlim()
    axs3.plot([start_pt, end_pt], [(start_pt*slope) + intercept, (end_pt*slope) + intercept], 'black')
    print('For overall fit for R375, slope = %.2f' % slope, ' and intercept = %.2f' % intercept)

    axs1.yaxis.set_label_coords(-0.3, 0.5)
    axs3.yaxis.set_label_coords(-0.1, 0.5)


    # statistics
    # res = stats.linregress(x, y)

    # Save and display
    figname = 'Fig9_fit_summary.pdf'
    plt.savefig(figname, format='pdf', pad_inches=0.1)
    plt.show()


if __name__ == '__main__':
    fig9_summary()
