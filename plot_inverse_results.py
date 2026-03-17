#  plot_inverse_results.py
#  This script takes the latest (from common_params.py) run and plots the results

import matplotlib.pyplot as plt
import scipy.stats as stats
from common_params import *  # import common values across all models
import subject_data


def plot_inverse_results(use_fwd_model, this_case, unsupervised,fit_mode):
    # Key constants to set
    save_fig = True
    plot_guess = False
    plot_ct_uncertainty = False
    ct_vals = []  # to get rid of a warning

    # Open file and load data
    if use_fwd_model:
        data_filename = INVOUTPUTDIR + this_case + '_fitResults_' + fit_mode+'.npy'
        [_, rposvals, survvals, thrsim, thrtargs, initvec, [fitrposvals, fitsurvvals],
         _, rpos_err_metric, survivalerrs] = np.load(data_filename, allow_pickle=True)
    else:
        ct_vals = subject_data.subj_ct_data(this_case)
        data_filename = INVOUTPUTDIR + this_case + '_fitResults_' + fit_mode+'.npy'
        [_, rposvals, survvals, thrsim, thrtargs, initvec, [fitrposvals, fitsurvvals],
         _, rpos_err_metric, _, ct_vals] = np.load(data_filename, allow_pickle=True)

    # Make plots
    xvals = np.arange(0, NELEC) + 1
    l_e = NELEC - 1  # last electrode to plot

    # All on one plot
    figrows = 3
    figcols = 1
    fig_consol, axs = plt.subplots(figrows, figcols)
    fig_consol.set_figheight(9)
    fig_consol.set_figwidth(7.5)

    if not tp_extend:
        axs[0].plot(xvals[1:l_e] + 0.1, thrsim[0][0][1:l_e], marker='o', color='lightblue', label='fit MP')
        axs[0].plot(xvals[1:l_e] - 0.1, thrtargs[0][0][1:l_e], marker='o', color='blue', label='measured MP')
        axs[0].plot(xvals[1:l_e] + 0.1, thrsim[1][0][1:l_e], marker='o', color='pink', label='fit TP')
        axs[0].plot(xvals[1:l_e] - 0.1, thrtargs[1][0][1:l_e], marker='o', color='red', label='measured TP')
        mean_thr_err = (np.nanmean(np.abs(np.array(thrsim[0]) - np.array(thrtargs[0]))) +
                        np.nanmean(np.abs(np.array(thrsim[1]) - np.array(thrtargs[1])))) / 2.0
    else:

        axs[0].plot(xvals + 0.1, thrsim[0][0], marker='o', color='lightblue', label='fit MP')
        axs[0].plot(xvals - 0.1, thrtargs[0][0], marker='o', color='blue', label='measured MP')
        axs[0].plot(xvals[1:l_e] + 0.1, thrsim[1][0][1:l_e], marker='o', color='pink', label='fit TP')
        axs[0].plot(xvals[1:l_e] - 0.1, thrtargs[1][0][1:l_e], marker='o', color='red', label='measured TP')
        mean_thr_err = (np.nanmean(np.abs(np.array(thrsim[0]) - np.array(thrtargs[0]))) +
                        np.nanmean(np.abs(np.array(thrsim[1]) - np.array(thrtargs[1])))) / 2.0

    yl = 'Threshold (dB)'
    if use_fwd_model:
        title_text = 'Known scenario thresholds: ' + this_case + '; mean thr error (dB): ' + '%.2f' % mean_thr_err
    else:
        title_text = ('Subject thresholds: ' + this_case + ' ' + RE_TEXT + ' ' + RI_TEXT + ' ' + TARG_TEXT +
                      '; mean thr error (dB): ' + '%.2f' % mean_thr_err)
    axs[0].set(xlabel='Electrode number', ylabel=yl, title=title_text)
    axs[0].set_xlim(0, 17)
    axs[0].legend(loc='upper right', ncol=2)

    if use_fwd_model:
        [dist_corr, dist_corr_p] = stats.pearsonr(1 - rposvals, 1 - fitrposvals)
    else:
        if np.any(ct_vals):
            [dist_corr, dist_corr_p] = stats.pearsonr(1 - ct_vals, 1 - fitrposvals)
    title_text = 'Fit and actual positions; mean position error (mm): ' + '%.2f' % rpos_err_metric
    if not tp_extend:

        axs[1].plot(xvals[1:l_e] + 0.1, 1 - fitrposvals[1:l_e], marker='o', color='gray', label='fit')
        if use_fwd_model:
            axs[1].plot(xvals[1:l_e] - 0.1, 1 - rposvals[1:l_e], marker='o', color='black', label='actual')
        else:
            if np.any(ct_vals):
                axs[1].plot(xvals[1:l_e] - 0.1, 1 - ct_vals[1:l_e], marker='o', color='black', label='CT estimate')
                if plot_ct_uncertainty:
                    axs[1].fill_between(xvals[1:l_e] - 0.1, (1 - ct_vals[1:l_e]) - ct_uncertainty, (1 - ct_vals[1:l_e]) + ct_uncertainty,
                                        color='black', alpha=0.1)
    else:
        axs[1].plot(xvals + 0.1, 1 - fitrposvals, marker='o', color='gray', label='fit')
        if use_fwd_model:
            axs[1].plot(xvals - 0.1, 1 - rposvals, marker='o', color='black', label='actual')
        else:
            if np.any(ct_vals):
                axs[1].plot(xvals - 0.1, 1 - ct_vals, marker='o', color='black', label='CT estimate')
                if plot_ct_uncertainty:
                    axs[1].fill_between(xvals - 0.1, (1 - ct_vals) - ct_uncertainty, (1 - ct_vals) + ct_uncertainty,
                                        color='black', alpha=0.1)

   # if use_fwd_model:
    #    axs[1].plot(xvals - 0.1, 1 - rposvals, marker='o', color='black', label='actual')
   # else:
    #    if np.any(ct_vals):
     #       axs[1].plot(xvals - 0.1, 1 - ct_vals, marker='o', color='black', label='CT estimate')
      #      if plot_ct_uncertainty:
       #         axs[1].fill_between(xvals - 0.1, (1 - ct_vals) - ct_uncertainty, (1 - ct_vals) + ct_uncertainty,
        #                            color='black', alpha=0.1)

    if plot_guess:
        axs[1].plot(xvals, 1 - initvec[0:NELEC], marker='o', color='purple', label='initial guess')

    axs[1].set(xlabel='Electrode number', ylabel='Electrode distance (mm)', title=title_text)
    axs[1].set_xlim(0, 17)
    axs[1].set_ylim(0, 2)
    axs[1].legend()

    # TODO conisider using a horizontal line rather than a dot to show that survival value is for a range of space.
    title_text = 'Fit survival values'
    if not tp_extend:
        if use_fwd_model:
            # axs[2].plot(xvals[1:l_e], fitsurvvals[1:l_e], marker='o', color='red', label='fit')
            for el in range(NELEC-2):
                axs[2].plot([xvals[el+1] - 0.45, xvals[el+1] + 0.45], [survvals[el+1], survvals[el+1]],
                            color='black', label='desired')
            axs[2].plot(xvals[1:l_e], survvals[1:l_e], marker='o', color='green', label='desired')
            if plot_guess:
                axs[2].plot(xvals[1:l_e], initvec[NELEC-2:], marker='o', color='purple', label='initial guess')

        else:
            # axs[2].plot(xvals[1:l_e], fitsurvvals[1:l_e], marker='o', color='black', label='modeled')
            for el in range(NELEC - 2):
                axs[2].plot([xvals[el+1] - 0.45, xvals[el+1] + 0.45], [fitsurvvals[el], fitsurvvals[el]],
                            color='black')

    else:

        if use_fwd_model:
            # axs[2].plot(xvals, fitsurvvals, marker='o', color='red', label='fit')
            for el in range(NELEC):
                axs[2].plot([xvals[el] - 0.45, xvals[el] + 0.45], [survvals[el], survvals[el]], marker='o', color='green', label='desired')
            if plot_guess:
                axs[2].plot(xvals, initvec[NELEC:], marker='o', color='purple', label='initial guess')

        else:
            axs[2].plot(xvals, fitsurvvals, marker='o', color='black', label='modeled')

    axs[2].set(xlabel='Electrode number', ylabel='Fractional neuronal density', title=title_text)
    axs[2].set_xlim(0, 17)
    axs[2].set_ylim(0, 1)
    axs[2].legend()
    fig_consol.tight_layout()

    # -- could add plots of error (difference between desired/measured and fitted values)

    if save_fig:
        save_file_name = INVOUTPUTDIR + this_case + '_fitResultsFig_' + fit_mode + '.png'
        fig_consol.savefig(save_file_name)

    # test correlation figure
    if not use_fwd_model:
        if np.any(ct_vals):
            fig2, ax2 = plt.subplots()
            ct_dist = 1 - ct_vals[1:l_e]
            fit_dist = 1 - fitrposvals[1:l_e]
            fig2 = plt.plot(ct_dist, fit_dist, 'o')
            ax2.set(xlabel='CT distance (mm)', ylabel='fit distance (mm)')
            ax2.set_xlim([0, 2])
            ax2.set_ylim([0, 2])
            fit_corr = np.corrcoef(ct_dist, fit_dist)
    if not unsupervised:
        plt.show()


    # second plot, scatter of electrode position and neuronal survival
    fig3, axs3 = plt.subplots(1, 1)
    axs3.plot(1 - fitrposvals, fitsurvvals, 'o', color='green')
    axs3.set_xlabel('fit distance (mm)')
    axs3.set_ylabel('fit density')
    axs3.set_xlim([0, 2])
    axs3.set_ylim([0, 1])
    # now plot initial guesses
    if not tp_extend:
        axs3.plot(1 - initvec[0:NELEC-2], initvec[NELEC-2:], 'o', color='red')
        for i in range(NELEC-2):
            axs3.plot([1 - fitrposvals[i], 1 - initvec[i]], [fitsurvvals[i], initvec[NELEC-2 + i]], color='black')

    else:
        axs3.plot(1 - initvec[0:NELEC], initvec[NELEC:], 'o', color='red')
        for i in range(NELEC):
            axs3.plot([1 - fitrposvals[i], 1 - initvec[i]], [fitsurvvals[i], initvec[NELEC + i]], color='black')

    title_text = this_case + '  Mean position error (mm): ' + '%.2f' % rpos_err_metric
    fig3.suptitle(title_text, fontsize=14)
    axs3.text(1.75, 0.95, 'guess', color='red')
    axs3.text(1.75, 0.9, 'fit', color='green')
    save_file_name = INVOUTPUTDIR + this_case + '_fitResults_scatter.png'
    fig3.savefig(save_file_name, format='png')
    if not unsupervised:
        plt.show()

if __name__ == '__main__':
    use_fwd_model = True
    case = scenarios[2]
    # case = 'S52'
    # espace = 0.85
    # case = 'RampRposSGradual80'
    # case = 'Gradual80R00'

    unsupervised = False
    fit_mode='combined'
    plot_inverse_results(use_fwd_model, case, unsupervised,fit_mode)
