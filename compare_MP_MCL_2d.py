# Script to compare 2D contours for monopolar threshold and MCL stimulations

# Import required packages
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from scipy import interpolate
from common_params import *


# noinspection PyUnusedLocal
def twod_mp_mcl_compare():

    ver = '2_5'
    # this is an option to search systematically for unique or multiple solutions for a given set of
    # monopolar and tripolar threshold values

    # filled_contours = False  # Fill the contour plots? Otherwise make contour lines
    # Set default figure values
    mpl.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2

    # # Declare variables  TOTO should get these from the saved or some common param
    surv_vals = np.arange(0.04, 0.97, 0.02)
    nsurv = len(surv_vals)
    rpos_vals = np.arange(-0.95, 0.96, 0.02)
    nrpos = len(rpos_vals)

    espace = 1.1
    sigval = 0.0
    thr_vals = [100, 500]
    if espace == 0.85:
        e_txt = '085'
    elif espace == 1.1:
        e_txt = '110'
    else:
        e_txt = 'xxx'
    es_text = '_espace_' + e_txt

    # Load first sigma data
    # Rename fwd and inverse output directories
    new_dir_suffix = 'RE%d' % R_EXT + '_'+'RI%d' % R_INT +'std_%.1f' % ACT_STDREL + '_thr_%d' % thr_vals[0] + '/'
    # offset = len(cp.FWD_OUT_PRFIX)
    FWD_OUT_PRFIX = 'FWD_OUTPUT/'
    FWDOUTPUTDIR = FWD_OUT_PRFIX + new_dir_suffix

    sigtext = f"{sigval:.1f}".replace(".", "_")

    datafile = FWDOUTPUTDIR + '2D_' + sigtext + '_' + STD_TEXT + es_text + '.csv'

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

    print("average monopolar threshold: ", np.nanmean(mono_thr))

    # Load mcl data
    new_dir_suffix = 'RE%d' % R_EXT + '_'+'RI%d' % R_INT +'std_%.1f' % ACT_STDREL + '_thr_%d' % thr_vals[1] + '/'
    FWD_OUT_PRFIX = 'FWD_OUTPUT/'
    FWDOUTPUTDIR = FWD_OUT_PRFIX + new_dir_suffix
    sigtext = f"{sigval:.1f}".replace(".", "_")
    datafile = FWDOUTPUTDIR + '2D_' + sigtext + '_' + STD_TEXT + es_text + '.csv'


    file = open(datafile)
    numlines = len(file.readlines())
    file.close()

    with open(datafile, newline='') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',')
        ncol = len(next(datareader))
        csvfile.seek(0)
        mcl_thr = np.empty([numlines, ncol])
        for i, row in enumerate(datareader):
            # Do the parsing
            mcl_thr[i, :] = row

    print('Data retrieved from directory: ', FWDOUTPUTDIR)

    # Measure min/max/mean differences between monopolar and tripolar
    thr_diff = mcl_thr - mono_thr
    mean_diff = np.mean(thr_diff[:])
    min_diff = np.min(thr_diff[:])
    max_diff = np.max(thr_diff[:])
    print('Min/max/mean differences: ', min_diff, ' , ', max_diff, ' , ', mean_diff)

    # # set up 2D interpolation
    rp_curt = rpos_vals[0:-2]
    xnew = np.linspace(rpos_vals[1], rpos_vals[-1], 50)
    ynew = np.linspace(surv_vals[1], surv_vals[-1], 50)
    np.meshgrid(xnew, ynew)

    f_interp = interpolate.RectBivariateSpline(rp_curt, surv_vals, mono_thr[:, 0:-2].transpose())
    znew_mp = f_interp(xnew, ynew)
    m_min = np.min(znew_mp)
    m_max = np.max(znew_mp)

    # f_interp = interpolate.interp2d(rp_curt, surv_vals, tripol_thr[:, 0:-2])
    f_interp = interpolate.RectBivariateSpline(rp_curt, surv_vals, mcl_thr[:, 0:-2].transpose())
    znew_tp = f_interp(xnew, ynew)
    t_min = np.min(znew_tp)
    t_max = np.max(znew_tp)

    # all_min = np.min([t_min, m_min])  # if we want to do this automatically
    # all_max = np.max([t_max, m_max])
    # rounding manually for figure
    all_min = 30.0
    all_max = 90.0
    n_levels = 6
    labels = ['P1', 'P2', 'P3', 'P4']

    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    fig.tight_layout(pad=3, w_pad=2, h_pad=2.0)
    plt.figtext(0.03, 0.93, 'A', fontsize=20, weight='bold')
    plt.figtext(0.51, 0.93, 'B', fontsize=20, weight='bold')


    cs3 = axs[0].contour(1 - rpos_vals, surv_vals, mono_thr,
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
    sigtext = 'MP thresh. Sigma: ' + f'{sigval:.1f}'
    axs[0].set_title(sigtext, fontsize=12)
    lab_shift = 0.025
    axs[0].text(1 - high_rpos_val - 0.15, high_surv_val - 0.06, labels[0], horizontalalignment='left',
                   verticalalignment='bottom', fontweight='bold')
    axs[0].text(1 - low_rpos_val + lab_shift, high_surv_val - 0.06, labels[1], horizontalalignment='left',
                   verticalalignment='bottom', fontweight='bold')
    axs[0].text(1 - high_rpos_val - 0.15, low_surv_val, labels[2], horizontalalignment='left',
                   verticalalignment='bottom', fontweight='bold')
    axs[0].text(1 - low_rpos_val + lab_shift, low_surv_val, labels[3], horizontalalignment='left',
                   verticalalignment='bottom', fontweight='bold')
    axs[0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [high_surv_val, high_surv_val], color='blue')
    axs[0].plot([np.min(1 - rpos_vals), np.max(1 - rpos_vals)], [low_surv_val, low_surv_val], color='blue',
                   linestyle='dashed')
    axs[0].plot([1 - low_rpos_val, 1 - low_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black',
                   linestyle='dashed')
    axs[0].plot([1 - high_rpos_val, 1 - high_rpos_val], [np.min(surv_vals), np.max(surv_vals)], color='black')
    axs[0].set_xticks([0.4, 0.8, 1.2, 1.6])


    # pos_1=np.array([0.476232742, 0.476232742, 0.575021851, 0.656158433, 0.745848506, 0.812965907, 0.856014243, 0.968904478, 1.047277507, 1.103322552, 1.19391566, 1.294313055, 1.366825349, 1.432090617, 1.51687377, 1.51687377])
    # surv_1=np.array([0.81395243, 0.81395243, 0.526326022, 0.429863864, 0.473558241, 0.553300381, 0.691186001, 0.986029868, 0.758011284, 0.434596852, 0.518933325, 0.477039127, 0.545469618, 0.721835387, 0.819748833, 0.819748833])
    # axs[0].plot(pos_1,surv_1, 'o')
    # all_min = 30.0
    # all_max = 90.0
    # n_levels = 6
    all_min = 25
    all_max = 75
    n_levels =10
    cs4 = axs[1].contour(1 - rpos_vals, surv_vals, mcl_thr,
                            np.arange(all_min, all_max, (all_max - all_min) / n_levels),
                            extend='both', colors='k')
    manual_locations = [(0, 0.9), (0.3, 0.3), (0.7, 0.7), (1.2, 0.2), (1.2, 0.5), (1.3, 0.02)]

    axs[1].clabel(cs4, cs4.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
    sigtext = 'MCL: 500 neurons. Sigma: ' + f'{sigval:.1f}'
    axs[1].set_title(sigtext, fontsize=12)
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
        #example values from uniform scenerio
       # pos_1=np.array([0.476232742, 0.476232742, 0.575021851, 0.656158433, 0.745848506, 0.812965907, 0.856014243, 0.968904478, 1.047277507, 1.103322552, 1.19391566, 1.294313055, 1.366825349, 1.432090617, 1.51687377, 1.51687377])
       # surv_1=np.array([0.81395243, 0.81395243, 0.526326022, 0.429863864, 0.473558241, 0.553300381, 0.691186001, 0.986029868, 0.758011284, 0.434596852, 0.518933325, 0.477039127, 0.545469618, 0.721835387, 0.819748833, 0.819748833])
       # axs[0, 1].plot(pos_1,surv_1, 'o')
#        pos_1=np.array([1.105082129, 1.105082129, 1.078369686, 0.894760897, 1.042783396, 1.050013736, 1.139597599, 1.226986237, 1.271727668, 1.200529353, 1.068024186, 0.7044402, 1.00113993, 1.153762216, 1.214020265, 1.214020265])
  #      surv_1=np.array([0.635876938, 0.635876938, 0.626996472, 0.245047452, 0.586458575, 0.676018414, 0.811287604, 0.830225343, 0.734071391, 0.998805076, 0.334925774, 0.000260672, 0.318277833, 0.758215672, 0.789004392, 0.789004392])
  #      axs[0, 1].plot(pos_1,surv_1, 'o')

    #  calculate  indices from target survival and rpos values
    low_rpos_idx = np.argmin(np.abs(rpos_vals - low_rpos_val))
    high_rpos_idx = np.argmin(np.abs(rpos_vals - high_rpos_val))
    low_surv_idx = np.argmin(np.abs(surv_vals - low_surv_val))
    high_surv_idx = np.argmin(np.abs(surv_vals - high_surv_val))


    filename = 'fig_twod_contour' + es_text + '.pdf'
    # plt.savefig('Fig4_2D_contour.eps', format='eps')
    plt.savefig(filename, format='pdf')

    fig2, axs2 = plt.subplots(1, 2, figsize=(12, 5))

    cs3 = axs2[0].contour(1 - rpos_vals, surv_vals, mono_thr,
                         np.arange(all_min, all_max, (all_max - all_min) / n_levels), colors='k')
    manual_locations = [(0, 0.9), (0.3, 0.3), (0.7, 0.7), (1.2, 0.2), (1.2, 0.5), (1.3, 0.02)]

    axs2[0].clabel(cs3, cs3.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)

    cs4 = axs2[0].contour(1 - rpos_vals, surv_vals, mcl_thr,
                     np.arange(all_min, all_max, (all_max - all_min) / n_levels), colors='r')

    manual_locations = [(0, 0.9), (0.3, 0.3), (0.7, 0.7), (1.2, 0.2), (1.2, 0.5), (1.3, 0.02)]

    axs2[0].clabel(cs4, cs4.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
    sigtext = ver + ' MP thresh. Sigma: ' + f'{sigval:.1f}'
    axs2[0].set_title(sigtext, fontsize=12)
    axs2[0].set_xlabel('Electrode distance (mm)')
    axs2[0].set_xticks([0.4, 0.8, 1.2, 1.6])
    axs2[0].set_ylabel('Fractional neuronal density')

    # Second panel with difference
    diff_vals = mcl_thr-mono_thr
    cs5 = axs2[1].contour(1 - rpos_vals, surv_vals, diff_vals, [2, 3, 4, 5, 6, 7, 8, 9.9, 10, 10.5, 11, 16])
    manual_locations = [(0.1, 0.1), (0.25, 0.6), (0.6, 0.2), (1.2, 0.2), (0.4, 0.4), (0.8, 0.3), (1.6, 0.6)]
    axs2[1].clabel(cs5, cs5.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
    sigtext = ver + ' MCL-MP thresh Difference. sigma =  ' + f'{sigval:.1f}'
    axs2[1].set_title(sigtext, fontsize=12)
    axs2[1].set_xlabel('Electrode distance (mm)')
    axs2[1].set_xticks([0.4, 0.8, 1.2, 1.6])
    axs2[1].set_ylabel('Fractional neuronal density')
    plt.show()



    plt.show()



if __name__ == '__main__':
    twod_mp_mcl_compare()
