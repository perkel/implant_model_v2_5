#  fig_inverse_results.py
#  This script takes the latest (from common_params.py) run and plots the results

import matplotlib.pyplot as plt
import matplotlib as mpl
from common_params import *  # import common values across all models


def fig_scenario_inverse_results():
    # Key constants to set
    figscenarios = ['Gradual80R00', 'RampRposS80', 'RampRposSGradual80']
    # Set default figure values
    mpl.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2

    # Directory names
    new_dir_suffix = 'R%d' % R_EXT + '_' + 'std_%.1f' % ACT_STDREL + '_thr_%d' % THRTARG
    # offset = len(cp.FWD_OUT_PRFIX)
    INV_OUT_PRFIX = 'INV_OUTPUT/'
    INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix


    # All on one plot
    figrows = 3
    figcols = 3
    ncols = 3
    fig_consol, axs = plt.subplots(figrows, figcols)
    fig_consol.set_figheight(8)
    fig_consol.set_figwidth(11)
    xvals = np.arange(0, NELEC) + 1
    l_e = NELEC - 1  # last electrode to plot
    plt.figtext(0.02, 0.96, 'A', fontsize=20, fontweight='bold')
    plt.figtext(0.35, 0.96, 'B', fontsize=20, fontweight='bold')
    plt.figtext(0.67, 0.96, 'C', fontsize=20, fontweight='bold')

    # Loop on datafiles
    for i, scenario in enumerate(figscenarios):

        ax = axs.flat[i]
        # Open file and load data
        data_filename = INVOUTPUTDIR + scenario + '_fitResults_' + 'combined.npy'

        data = np.load(data_filename, allow_pickle=True)
        [_, rposvals, survvals, thrsim, thrtargs, _, [fitrposvals, fitsurvvals], _, _, _] = data

        # Make plots
        ax.plot(xvals[1:l_e] - 0.2, thrsim[1][0][1:l_e], marker='^', color='purple', label='TP fit')
        ax.plot(xvals[1:l_e] + 0.2, thrtargs[1][0][1:l_e], marker='^', color='orange', label='TP actual')
        ax.plot(xvals - 0.2, thrsim[0][0], marker='o', color='purple', label='MP fit')
        ax.plot(xvals + 0.2, thrtargs[0][0], marker='o', color='orange', label='MP actual')
        # title_text = scenario
        # ax.set_title(title_text)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        xlim = ax.get_xlim()
        ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
        ax.set_ylim(50, 92)
        if i == 0:
            ax.set_ylabel('Threshold (dB re 1 $\mu$A)', font='Arial')

        if i == 1:
            ax.legend(framealpha=0.0, fontsize=10)

        ax = axs.flat[i + ncols]
        # Note converting rpos to distance from inner wall, calling it "Electrode distance"
        ax.plot(xvals[1:l_e] + 0.2, 1-rposvals[1:l_e], marker='o', color='orange', label='actual')
        ax.plot(xvals[1:l_e] - 0.2, 1-fitrposvals[1:l_e], marker='o', color='purple', label='estimated')

        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xlim(xlim)
        ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
        ax.set_ylim(0.2, 1.7)
        if i == 0:
            ax.set(ylabel='Electrode distance (mm)')
            ax.legend(framealpha=0.0, fontsize=10)

        ax = axs.flat[i + 2*ncols]
        ax.plot(xvals[1:l_e] + 0.2, survvals[1:l_e], marker='o', color='orange', label='actual')
        ax.plot(xvals[1:l_e] - 0.2, fitsurvvals[1:l_e], marker='o', color='purple', label='estimated')
        ax.set(xlabel='Electrode number')
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xlim(xlim)
        ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16])
        ax.set_ylim(0.25, 0.95)
        if i == 0:
            ax.set(ylabel='Fractional neuronal density')
            # ax.legend()
        if i == ncols - 1:
            ax.set(xlabel='Electrode number')

        fig_consol.tight_layout()

        # -- could add plots of error (difference between desired/measured and fitted values)
        # Calculate error values
        max_mp_error = np.max(np.abs(np.subtract(thrsim[0][0], thrtargs[0][0])))
        max_tp_error = np.max(np.abs(np.subtract(thrsim[1][0][1:l_e], thrtargs[1][0][1:l_e])))
        max_dist_error = np.max(np.abs(np.subtract(rposvals[1:l_e], fitrposvals[1:l_e])))
        max_density_error = np.max(np.abs(np.subtract(survvals, fitsurvvals)))
        print("scenario: ", scenario, " max mp, tp and dist and surv errors: ",
              max_mp_error, max_tp_error, max_dist_error, max_density_error)

    figname = "fig6_scenarios.pdf"
    fig_consol.savefig(figname, format='pdf')
    figname = "fig6_scenarios.png"
    fig_consol.savefig(figname, format='png')
    plt.show()


if __name__ == '__main__':
    fig_scenario_inverse_results()
