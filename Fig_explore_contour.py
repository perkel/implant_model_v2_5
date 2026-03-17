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


def read_inv_summary(res, scen):
    # Reads a summary file, and tests whether average rpos error is less than chance based on shuffling
    # You need to run the inverse model for all subjects before making this figure

    # construct correct path for this resistivity
    R_TEXT = 'R' + str(round(res))
    INV_OUT_PRFIX = 'INV_OUTPUT/'
    INVOUTPUTDIR = INV_OUT_PRFIX + R_TEXT + ACTR_TEXT + STD_TEXT + TARG_TEXT

    summary_file_name = INVOUTPUTDIR + scen + '_fitResults_combined.npy'
    data = np.load(summary_file_name, allow_pickle=True)


    return

def fig_explore_contour():
    # Constants
    label_ypos = 1.05
    n_subj = 18
    nscen = len(scenarios)
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
    r_val = 375.0

    # Layout figure
    fig1, axs1 = plt.subplots(3, 2, figsize=(8, 8))
    fig1.tight_layout(pad=3)

    plt.figtext(0.02, 0.95, 'A', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.95, 'B', color='black', size=20, weight='bold')
    plt.figtext(0.02, 0.63, 'C', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.63, 'D', color='black', size=20, weight='bold')
    plt.figtext(0.02, 0.33, 'E', color='black', size=20, weight='bold')
    plt.figtext(0.49, 0.33, 'F', color='black', size=20, weight='bold')

    # get data
    # thr_sum_all_1, thresh_err_summary_1, rpos_fit_vals_1, rpos_vals_1, rpos_err_summary_1, aaa_temp, dist_corr_1, dist_corr_p_1 =(
    #     read_inv_summary(r_vals[1]))



    color = iter(cm.rainbow(np.linspace(0, 1, n_subj)))
    for idx, scen in enumerate(scenarios):  # Panel A
        read_inv_summary(r_val, scen)


        axs1[0, 0].plot(x, y, '.', color=fig9_colors[idx, :])
        axs1[0, 0].set_xlabel('Measured monopolar threshold (dB)')
        axs1[0, 0].set_ylabel('Fit monopolar threshold (dB)')
        axs1[0, 0].spines['top'].set_visible(False)
        axs1[0, 0].spines['right'].set_visible(False)
        axs1[0, 0].set_xlim([70, 90])
        axs1[0, 0].set_ylim([70, 90])



    # Save and display
    figname = 'Fig9_fit_summary.pdf'
    plt.savefig(figname, format='pdf', pad_inches=0.1)
    plt.show()


if __name__ == '__main__':
    fig_explore_contour()
