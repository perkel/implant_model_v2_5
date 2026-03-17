# Script to make figures 4 (3D plot) & 5 (example contours
# for Perkel, Goldwyn & Arenberg (2023) implant model paper on forward and inverse models

# Import required packages
import matplotlib.pyplot as plt
from scipy import interpolate
from common_params import *


def ecap_pars(stim, amp, start_amp, end_amp):
    min_x = np.argmax(amp > start_amp)
    max_x = np.argmax(amp > end_amp)

    temp_x = stim[min_x:max_x]
    # temp1 = np.array(ecaps[elecidx][segment][:
    temp_y = amp[min_x:max_x]

    #  adjust start and end points to reduce jitter
    my_x0 = np.interp(start_amp, [amp[min_x - 1], temp_y[0]], [stim[min_x - 1], temp_x[0]])
    my_x1 = np.interp(end_amp, [temp_y[-1], amp[max_x + 1]], [temp_x[-1], stim[max_x + 1]])
    temp_x2 = np.append(temp_x, my_x1)
    temp_x3 = np.append(my_x0, temp_x2)
    temp_y2 = np.append(temp_y, end_amp)
    temp_y3 = np.append(start_amp, temp_y2)
    # print('st & end amp: ', start_amp, end_amp, ' min_x, max_x: ', min_x, max_x, ' tx, ty: ',
    # '%.2f' % temp_x[0], '%.2f' % temp_y[0])
    slope, y_intercept = np.polyfit(temp_x3, temp_y3, 1)
    x_intercept = -y_intercept / slope
    return slope, x_intercept, y_intercept


# noinspection PyUnusedLocal


def fig2d_ecap_contour():
    # this is an option to search systematically for ecap thresholds and slopes
    # across electrode distance and neuronal density
    map_unique_solutions = False
    filled_contours = True  # Fill the contour plots? Otherwise, make contour lines
    # Set default figure values
    plt.style.use('matplotlib_settings')

    # Declare variables -- should be able to get these from the data files
    surv_vals = np.arange(0.04, 0.97, 0.02)
    nsurv = len(surv_vals)
    rpos_vals = np.arange(-0.95, 0.96, 0.02)
    nrpos = len(rpos_vals)

    espace = 1.1
    if espace == 0.85:
        e_txt = '085'
    elif espace == 1.1:
        e_txt = '110'
    else:
        e_txt = 'xxx'
    es_text = '_espace_' + e_txt

    # Load ecap data from binary file
    # Rename fwd and inverse output directories
    new_dir_suffix = 'R%d' % R_EXT + '_' + 'std_%.1f' % ACT_STDREL
    # offset = len(cp.FWD_OUT_PRFIX)
    FWD_OUT_PRFIX = 'FWD_OUTPUT/'
    # FWDOUTPUTDIR = FWD_OUT_PRFIX + new_dir_suffix

    # FWDOUTPUTDIR + 'ecap_' + STD_TEXT + es_text

    datafile = FWDOUTPUTDIR + "ecap_" + STD_TEXT + es_text + ".npz"
    ecap_data = np.load(datafile)
    surv_vals = ecap_data['arr_0']
    rpos_vals = ecap_data['arr_1']
    ecap_all = ecap_data['arr_2']

    l_marg = 10
    r_marg = 5

    ec2 = ecap_all[l_marg:-r_marg, l_marg:-r_marg, :, :]

    nsurv = len(surv_vals) - (l_marg + r_marg)
    nrpos = len(rpos_vals) - (l_marg + r_marg)
    slopes = np.zeros((nsurv, nrpos))
    intercepts = np.zeros((nsurv, nrpos))

    for sidx, s in enumerate(surv_vals[l_marg:-r_marg]):
        for ridx, r in enumerate(rpos_vals[l_marg:-r_marg]):
            # print('sidx, ridx:', sidx, ridx)
            stims = ec2[sidx, ridx, :, 0]
            amps = ec2[sidx, ridx, :, 1]
            end_amp = np.max(amps)
            st_amp = end_amp*0.9
            # st_amp = np.max([130, np.min(amps)])  # start amplitude

            ts, tix, tiy = ecap_pars(stims, amps, st_amp, end_amp)
            # plt.plot(stims, amps, '.')
            # x_vals = [tix, 60]
            # y_vals = np.zeros(2)
            # y_vals[0] = (ts*x_vals[0]) + tiy
            # y_vals[1] = (ts*x_vals[1]) + tiy
            # plt.plot(x_vals, y_vals, '-r')
            # plt.show()
            slopes[sidx, ridx] = ts
            intercepts[sidx, ridx] = tix

    print('Data retrieved from directory: ', FWDOUTPUTDIR)

    # fig0, a0 = plt.subplots(2, 1, sharex=True)
    # a0[0].imshow(slopes)
    # a0[1].imshow(intercepts)
    # plt.show()
    # # set up 2D interpolation
    rp_curt = rpos_vals[l_marg:-r_marg]
    sv_curt = surv_vals[l_marg:-r_marg]
    xnew = np.linspace(rpos_vals[1], rpos_vals[-1], 500)
    ynew = np.linspace(surv_vals[1], surv_vals[-1], 500)
    np.meshgrid(xnew, ynew)

    f_interp = interpolate.RectBivariateSpline(rp_curt, sv_curt, slopes.transpose())
    znew_sl = f_interp(xnew, ynew)
    sl_min = np.min(znew_sl)
    sl_max = np.max(znew_sl)

    # f_interp = interpolate.interp2d(rp_curt, surv_vals, tripol_thr[:, 0:-2])
    f_interp = interpolate.RectBivariateSpline(rp_curt, sv_curt, intercepts.transpose())
    znew_int = f_interp(xnew, ynew)
    int_min = np.min(znew_int)
    int_max = np.max(znew_int)

    # all_min = np.min([t_min, m_min])  # if we want to do this automatically
    # all_max = np.max([t_max, m_max])
    # rounding manually for figure
    # all_min = 10.0
    # all_max = 50.0
    all_min = np.floor(np.min(slopes[:]))
    all_max = np.ceil(np.max(slopes[:]))

    n_levels = 6
    labels = ['P1', 'P2', 'P3', 'P4']
    low_rpos_val = -0.5
    high_rpos_val = 0.5
    low_surv_val = 0.4
    high_surv_val = 0.8

    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    fig.tight_layout(pad=3, w_pad=2, h_pad=2.0)
    # plt.figtext(0.03, 0.93, 'A', fontsize=20, weight='bold')
    # plt.figtext(0.51, 0.93, 'B', fontsize=20, weight='bold')
    # plt.figtext(0.03, 0.49, 'C', fontsize=20, weight='bold')
    # plt.figtext(0.51, 0.49, 'D', fontsize=20, weight='bold')

    if filled_contours:
        # print('test')
        cs3 = axs[0].contourf(1 - rpos_vals[l_marg:-r_marg], surv_vals[l_marg:-r_marg], slopes,
                              np.arange(all_min, all_max, (all_max - all_min) / n_levels),
                              cmap='viridis', extend='both')
        axs[0].set_xlabel('Electrode distance (mm)')
        axs[0].set_ylabel('Fractional neuronal density')
        axs[0].set_title('Slopes', fontsize=12)
        cb3 = fig.colorbar(cs3, extend='neither', extendfrac=None)

        # low_rpos_val = -0.5  # these are the P1-4 points
        # high_rpos_val = 0.5
        # low_surv_val = 0.4
        # high_surv_val = 0.8
        # axs[0, 0].set_xlabel('Electrode distance (mm)')
        # axs[0, 0].set_ylabel('Fractional neuronal density')
        # axs[0, 0].set_title('Monopolar', fontsize=12)
        # lab_shift = 0.025
        # axs[0, 0].text(1 - high_rpos_val + lab_shift, high_surv_val, labels[0], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 0].text(1 - low_rpos_val + lab_shift, high_surv_val, labels[1], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 0].text(1 - high_rpos_val + lab_shift, low_surv_val, labels[2], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 0].text(1 - low_rpos_val + lab_shift, low_surv_val, labels[3], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 0].plot([1 - high_rpos_val], [high_surv_val], 'sk', markersize=20)
        # axs[0, 0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [high_surv_val, high_surv_val], color='blue')
        # axs[0, 0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [low_surv_val, low_surv_val], color='blue',
        #                linestyle='dashed')
        # axs[0, 0].plot([1 - low_rpos_val, 1 - low_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black',
        #                linestyle='dashed')
        # axs[0, 0].plot([1 - high_rpos_val, 1 - high_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black')
        # axs[0, 0].set_xticks([0.4, 0.8, 1.5, 1.6])
        #
        all_min = np.floor(np.min(intercepts[:]))
        all_max = np.ceil(np.max(intercepts[:]))
        cs4 = axs[1].contourf(1 - rpos_vals[l_marg:-r_marg], surv_vals[l_marg:-r_marg], intercepts,
                              np.arange(all_min, all_max, (all_max - all_min) / n_levels),
                              cmap='viridis', extend='both')
        axs[1].set_xlabel('Electrode distance (mm)')
        # axs[1].set_ylabel('Fractional neuronal density')
        axs[1].set_title('Thresholds (x-intercepts)', fontsize=12)
        cb4 = fig.colorbar(cs4, extend='neither', extendfrac=None)

        # axs[0, 1].set_title('Tripolar', fontsize=12)
        # axs[0, 1].set_xlabel('Electrode distance (mm)')
        # axs[0, 1].text(1 - high_rpos_val + lab_shift, high_surv_val, labels[0], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 1].text(1 - low_rpos_val + lab_shift, high_surv_val, labels[1], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 1].text(1 - high_rpos_val + lab_shift, low_surv_val, labels[2], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 1].text(1 - low_rpos_val + lab_shift, low_surv_val, labels[3], horizontalalignment='left',
        #                verticalalignment='bottom')
        # axs[0, 1].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [high_surv_val, high_surv_val], color='red')
        # axs[0, 1].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [low_surv_val, low_surv_val], color='red',
        #                linestyle='dashed')
        # axs[0, 1].plot([1 - low_rpos_val, 1 - low_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='gray',
        #                linestyle='dashed')
        # axs[0, 1].plot([1 - high_rpos_val, 1 - high_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='gray')
        # axs[0, 1].set_xticks([0.4, 0.8, 1.2, 1.6])
        plt.show()
    else:  # now for unfilled contours
        cs3 = axs[0].contour(1 - rpos_vals[l_marg:-r_marg], surv_vals[l_marg:-r_marg], slopes,
                             np.arange(all_min, all_max, (all_max - all_min) / n_levels), colors='k')

        # This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
        # then adds a percent sign.
        def fmt(x):
            s = f"{x:.1f}"
            if s.endswith("0"):
                s = f"{x:.0f}"
            return rf"{s} dB" if plt.rcParams["text.usetex"] else f"{s} dB"

        manual_locations = [(0, 0.9), (0.3, 0.6), (0.8, 0.3), (1.0, 0.2)]
        axs[0].clabel(cs3, cs3.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
        low_rpos_val = -0.5
        high_rpos_val = 0.5
        low_surv_val = 0.4
        high_surv_val = 0.8
        axs[0].set_xlabel('Electrode distance (mm)')
        axs[0].set_ylabel('Fractional neuronal density')
        axs[0].set_title('Slopes', fontsize=12)
        lab_shift = 0.025
        # axs[0, 0].text(1 - high_rpos_val - 0.15, high_surv_val - 0.06, labels[0], horizontalalignment='left',
        #                verticalalignment='bottom', fontweight='bold')
        # axs[0, 0].text(1 - low_rpos_val + lab_shift, high_surv_val - 0.06, labels[1], horizontalalignment='left',
        #                verticalalignment='bottom', fontweight='bold')
        # axs[0, 0].text(1 - high_rpos_val - 0.15, low_surv_val, labels[2], horizontalalignment='left',
        #                verticalalignment='bottom', fontweight='bold')
        # axs[0, 0].text(1 - low_rpos_val + lab_shift, low_surv_val, labels[3], horizontalalignment='left',
        #                verticalalignment='bottom', fontweight='bold')
        # axs[0, 0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [high_surv_val, high_surv_val], color='blue')
        # axs[0, 0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [low_surv_val, low_surv_val], color='blue',
        #                linestyle='dashed')
        # axs[0, 0].plot([1 - low_rpos_val, 1 - low_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black',
        #                linestyle='dashed')
        # axs[0, 0].plot([1 - high_rpos_val, 1 - high_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black')
        axs[0].set_xticks([0.4, 0.8, 1.2, 1.6])

        # all_min = 30.0
        # all_max = 90.0
        all_min = np.min(intercepts[:])
        all_max = np.max(intercepts[:])

        n_levels = 6
        cs4 = axs[1].contour(1 - rpos_vals[10:-5], surv_vals[10:-5], intercepts,
                             np.arange(all_min, all_max, (all_max - all_min) / n_levels),
                             extend='both', colors='k')
        manual_locations = [(0, 0.9), (0.3, 0.3), (0.7, 0.7), (1.2, 0.2), (1.2, 0.5), (1.3, 0.02)]

        # axs[1].clabel(cs4, cs4.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
        axs[1].set_title('Intercept', fontsize=12)
        axs[1].set_xlabel('Electrode distance (mm)')
        axs[1].text(1 - high_rpos_val + 0.03, high_surv_val - 0.06, labels[0], horizontalalignment='left',
                    verticalalignment='bottom', fontweight='bold')
        axs[1].text(1 - low_rpos_val + lab_shift, high_surv_val - 0.06, labels[1], horizontalalignment='left',
                    verticalalignment='bottom', fontweight='bold')
        axs[1].text(1 - high_rpos_val + 0.06, low_surv_val, labels[2], horizontalalignment='left',
                    verticalalignment='bottom', fontweight='bold')
        axs[1].text(1 - low_rpos_val + lab_shift, low_surv_val, labels[3], horizontalalignment='left',
                    verticalalignment='bottom', fontweight='bold')
        axs[1].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [high_surv_val, high_surv_val], color='red')
        axs[1].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [low_surv_val, low_surv_val], color='red',
                    linestyle='dashed')
        axs[1].plot([1 - low_rpos_val, 1 - low_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='gray',
                    linestyle='dashed')
        axs[1].plot([1 - high_rpos_val, 1 - high_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='gray')
        axs[1].set_xticks([0.4, 0.8, 1.2, 1.6])
        plt.show()

    # End contours

    #  calculate  indices from target survival and rpos values
    low_rpos_idx = np.argmin(np.abs(rpos_vals - low_rpos_val))
    high_rpos_idx = np.argmin(np.abs(rpos_vals - high_rpos_val))
    low_surv_idx = np.argmin(np.abs(surv_vals - low_surv_val))
    high_surv_idx = np.argmin(np.abs(surv_vals - high_surv_val))

    # axs[1, 0].plot(1 - rpos_vals, mono_thr[high_surv_idx, :], color='blue', linestyle='solid')
    # axs[1, 0].plot(1 - rpos_vals, tripol_thr[high_surv_idx, :], color='red', linestyle='solid')
    # axs[1, 0].plot(1 - rpos_vals, mono_thr[low_surv_idx, :], color='blue', linestyle='dashed')
    # axs[1, 0].plot(1 - rpos_vals, tripol_thr[low_surv_idx, :], color='red', linestyle='dashed')
    # axs[1, 0].axes.set_xlabel('Electrode distance (mm)')
    # axs[1, 0].axes.set_ylabel('Threshold (dB)')
    # axs[1, 0].axes.set_xlim([0.1, 1.9])
    # axs[1, 0].axes.set_ylim([30, 80])
    # axs[1, 0].set_xticks([0.4, 0.8, 1.2, 1.6])
    #
    # axs[1, 1].plot(surv_vals, mono_thr[:, high_rpos_idx], color='black', linestyle='solid')
    # axs[1, 1].plot(surv_vals, tripol_thr[:, high_rpos_idx], color='gray', linestyle='solid')
    # axs[1, 1].plot(surv_vals, mono_thr[:, low_rpos_idx], color='black', linestyle='dashed')
    # axs[1, 1].plot(surv_vals, tripol_thr[:, low_rpos_idx], color='gray', linestyle='dashed')
    # axs[1, 1].axes.set_xlabel('Fractional neuronal density')
    # axs[1, 1].axes.set_xlim([0.1, 0.9])
    # axs[1, 1].axes.set_ylim([30, 80])

    # filename = 'Fig4_2D_contour' + es_text + '.pdf'
    # # plt.savefig('Fig4_2D_contour.eps', format='eps')
    # plt.savefig(filename, format='pdf')

    figrows = 3
    figcols = 3
    fig2, axs2 = plt.subplots(figrows, figcols)

    for row in range(figrows):
        for col in range(figcols):
            eshape = ec2.shape
            s_idx = int(np.floor((row + 1) * eshape[0] / figrows) - 1)
            r_idx = int(np.floor(col + 1) * eshape[1] / figcols) - 1
            print('s_idx, r_idx: ', s_idx, r_idx)
            stims = ec2[s_idx, r_idx, :, 0]
            amps = ec2[s_idx, r_idx, :, 1]
            axs2[row, col].plot(stims, amps, '.r')
    plt.show()

    # idx_surv = 0
    # idx_rpos = 0
    # the_ax = 0
    # for i in range(4):
    #     if i == 0:
    #         idx_surv = high_surv_idx
    #         idx_rpos = high_rpos_idx
    #         the_ax = axs2[0, 0]
    #     elif i == 1:
    #         idx_surv = high_surv_idx
    #         idx_rpos = low_rpos_idx
    #         the_ax = axs2[0, 1]
    #     elif i == 2:
    #         idx_surv = low_surv_idx
    #         idx_rpos = high_rpos_idx
    #         the_ax = axs2[1, 0]
    #     elif i == 3:
    #         idx_surv = low_surv_idx
    #         idx_rpos = low_rpos_idx
    #         the_ax = axs2[1, 1]
    #
    #     this_mp_thr = [mono_thr[idx_surv, idx_rpos]]
    #     print('this_mp_thr: ', i, ' and ', this_mp_thr)
    #     cont_mp = the_ax.contour(1 - rpos_vals, surv_vals, mono_thr, this_mp_thr, colors='blue')
    #     if i == 2 or i == 3:
    #         the_ax.axes.set_xlabel('Electrode distance (mm)', fontsize=14)
    #
    #     if i == 0:
    #         the_ax.axes.set_ylabel('Fractional neuronal density', fontsize=14)
    #         the_ax.yaxis.set_label_coords(-0.2, -0.08)
    #
    #     this_tp_thr = [tripol_thr[idx_surv, idx_rpos]]
    #     print('this_tp_thr: ', i, ' and ', this_tp_thr)
    #     cont_tp = the_ax.contour(1 - rpos_vals, surv_vals, tripol_thr, this_tp_thr, colors='red')
    #     mpcontour = cont_mp.allsegs[0]
    #     tpcontour = cont_tp.allsegs[0]
    #     nmp = len(mpcontour[0])
    #     ntp = len(tpcontour[0])
    #     mpx = np.zeros(nmp)
    #     mpy = np.zeros(nmp)
    #     tpx = np.zeros(ntp)
    #     tpy = np.zeros(ntp)
    #
    #     for j in range(0, nmp):  # Should be able to do this without for loops
    #         mpx[j] = mpcontour[0][j][0]
    #         mpy[j] = mpcontour[0][j][1]
    #
    #     for j in range(0, ntp):
    #         tpx[j] = tpcontour[0][j][0]
    #         tpy[j] = tpcontour[0][j][1]
    #
    #     x, y = intsec.intersection(mpx, mpy, tpx, tpy)  # find intersection(s)
    #     if i > 0 and len(x) > 0:
    #         the_ax.plot(x[-1], y[-1], 'x', color='black', markersize='10', mew=2.5)
    #     the_ax.plot(x, y, 'x', color='black', markersize='10', mew=2.5)
    #     the_ax.set_xlim([0, 1.9])
    #     the_ax.text(0.1, 0.8, labels[i], fontsize=16)
    #
    # plt.savefig('Fig5_contour_examples.pdf', format='pdf')
    #
    # if map_unique_solutions:
    #     n_sols = np.zeros((nrpos, nsurv), dtype=int)
    #     # Map whether there is only a single solutions giving these values (or 0 or > 1)
    #     for sidx, surv in enumerate(surv_vals):
    #         # print('sidx = ', sidx)
    #         print('Approximately ', 100 * (sidx * nrpos) / (nrpos * nsurv), ' % done.')
    #         for ridx, rpos in enumerate(rpos_vals):
    #             # print('ridx is ', ridx)
    #
    #             fig, ax1 = plt.subplots()
    #             ax1 = plt.contour(rpos_vals, surv_vals, mono_thr, [mono_thr[sidx, ridx]], colors='green')
    #             ax2 = plt.contour(rpos_vals, surv_vals, tripol_thr, [tripol_thr[sidx, ridx]], colors='red')
    #             mpcontour = ax1.allsegs[0]
    #             tpcontour = ax2.allsegs[0]
    #             # if ifPlotContours == False:
    #             #    plt.close(fig)
    #
    #             nmp = len(mpcontour[0])
    #             ntp = len(tpcontour[0])
    #             mpx = np.zeros(nmp)
    #             mpy = np.zeros(nmp)
    #             tpx = np.zeros(ntp)
    #             tpy = np.zeros(ntp)
    #
    #             for j in range(0, nmp):  # Should be able to do this without for loops
    #                 mpx[j] = mpcontour[0][j]
    #                 mpy[j] = mpcontour[0][j][1]
    #
    #             for j in range(0, ntp):
    #                 tpx[j] = tpcontour[0][j][0]
    #                 tpy[j] = tpcontour[0][j][1]
    #
    #             # try finding intersection using Shapely class
    #             if len(mpcontour[0]) == 1 or len(tpcontour[0]) == 1:
    #                 n_sols[ridx, sidx] = 0
    #             else:
    #                 ('len contours: ', len(mpcontour[0]), len(tpcontour[0]))
    #                 line1 = shap.LineString(mpcontour[0])
    #                 line2 = shap.LineString(tpcontour[0])
    #                 x = shap.intersection(line1, line2)
    #
    #                 # How many intersections? 0, 1 or more?  If single intersection use those values
    #                 if x.geom_type == 'Point':  # test for multiple identical values
    #
    #                     n_sols[ridx, sidx] = 1
    #
    #                 elif x.geom_type == 'MultiPoint':
    #                     print('Multipoint at sidx = ', sidx, ' and ridx: ', ridx, ' : ', x)
    #                     n_sols[ridx, sidx] = 2
    #                 else:
    #                     n_sols[ridx, sidx] = 0
    #
    #     solmap_file = FWDOUTPUTDIR + 'solution_map' + STD_TEXT
    #     np.savez(solmap_file, surv_vals, rpos_vals, n_sols)
    #     fig_nsols, ax_nsols = plt.subplots()
    #     ax_nsols.contour(rpos_vals, surv_vals, n_sols, [1, 2])
    #
    # figname = 'Fig5_contour_examples.pdf'
    # plt.savefig(figname, format='pdf', pad_inches=0.1)
    # figname = 'Fig5_contour_examples.png'
    # plt.savefig(figname, format='png', pad_inches=0.1)

    plt.show()


if __name__ == '__main__':
    fig2d_ecap_contour()
