#  Script to help check for any changes in results after "cosmetic" or "cleanup" work

import pandas as pd
from common_params import *
stringency = 0.0001


def is_scenario(sc):  # test whether this is a scenario or subject
    # if this scenario is a subject, set use_forward_model to be false
    if (sc[0] == 'A' or sc[0] == 'S') and sc[1:3].isnumeric():
        return False  # it's a subject
    else:
        return True


def max_frac_diff(f_ref, f_test):

    # Load values -- customized depending on type
    data_ref = pd.read_csv(f_ref)
    data_test = pd.read_csv(f_test)

    #  Calculate maximum % difference in values
    return np.max(np.abs(np.divide(np.subtract(data_test, data_ref), data_ref)))


# set up list of files for forward model
def compare_inv_npy(f_ref, f_test):

    # Load values -- customized depending on type
    d_ref = np.load(f_ref, allow_pickle=True)
    d_test = np.load(f_test, allow_pickle=True)
    #  Calculate maximum % difference in values
    if np.max(d_ref[3]) == 0.0:  # Thresholds
        max_t_d = np.max(np.abs(np.subtract(d_test[3], d_ref[3])))
        print('denom is zero, returning absolute difference')
    else:
        max_t_d = np.max(np.abs(np.divide(np.subtract(d_test[3], d_ref[3]), d_ref[3])))

    # print('Positions....')
    if np.max(d_ref[6][0]) == 0.0:  # Thresholds
        max_r_d = np.max(np.abs(np.subtract(d_test[6], d_ref[6])))
        print('denom is zero, returning absolute difference')
    else:
        max_r_d = np.max(np.abs(np.divide(np.subtract(d_test[6], d_ref[6]), d_ref[6])))

    return max_t_d, max_r_d


# Main script to check output between a reference set of results and a new set after code modifications
ref_dir = 'MainPaperData/'

# Compare forward output results (for scanarios not subjects)
for scen in scenarios:
    if is_scenario(scen):
        file = 'FwdModelOutput_' + scen + '.csv'
        ref_file = FWD_OUT_PRFIX + ref_dir + new_dir_suffix + file
        test_file = FWDOUTPUTDIR + file

        max_diff = max_frac_diff(ref_file, test_file)
        if max_diff > stringency:
            print('ALERT! Detected a percentage difference of: ', max_diff*100, '%')
            # then indicate where the differences are
        else:
            print('Excellent! No differences above threshold were found for scenario: ', scen)

# Now compare output of 2D forward model
for espace in [0.85, 1.1]:
    if espace == 0.85:
        e_txt = '085'
    elif espace == 1.1:
        e_txt = '110'
    else:
        e_txt = 'xxx'
    es_text = '_espace_' + e_txt

    # set up filename
    min_surv = 0.04
    max_surv = 0.96
    min_rpos = -0.95
    max_rpos = 0.95
    hires = '_hi_res'

    descrip = 'surv_%.2f' % min_surv + '_%.2f' % max_surv + "_rpos_%.2f" %\
              min_rpos + '_%.2f' % max_rpos + hires

    # first monopolar
    file = 'Monopolar_2D_' + STD_TEXT + es_text + '.csv'
    ref = 'MainPaperData/'
    ref_file = FWD_OUT_PRFIX + ref + new_dir_suffix + file
    test_file = FWDOUTPUTDIR + file

    print('Checking 2D forward model...')
    max_diff = max_frac_diff(ref_file, test_file)
    if max_diff > stringency:
        print('ALERT! Detected a percentage difference of: ', max_diff * 100, '%')
        # then indicate where the differences are
    else:
        print('Excellent! No differences above threshold were found for 2D monopolar.')

    # next tripolar
    file = 'Tripolar_09_2D_' + STD_TEXT + es_text + '.csv'
    ref = 'MainPaperData/'
    ref_file = FWD_OUT_PRFIX + ref + new_dir_suffix + file
    test_file = FWDOUTPUTDIR + file

    max_diff = max_frac_diff(ref_file, test_file)
    if max_diff > stringency:
        print('ALERT! Detected a percentage difference of: ', max_diff * 100, '%')
        # then indicate where the differences are
    else:
        print('Excellent! No differences above threshold were found for 2D Triploar.')

# set up list of files for inverse model
print('Checking inverse results...')
for i, scen in enumerate(scenarios):
    print('Examining results for scenario/subject: ', scen)
    file = scen + '_fitResults_combined.npy'  # First do npy files
    ref = 'MainPaperData/'
    ref_file = INV_OUT_PRFIX + ref + new_dir_suffix + file
    test_file = INVOUTPUTDIR + file

    max_thresh_diff, max_rpos_diff = compare_inv_npy(ref_file, test_file)
    if max_thresh_diff > stringency:
        print('ALERT! Detected a threshold percentage difference of: ', max_thresh_diff*100, '%')
        # then indicate where the differences are
    else:
        print('Excellent! No differences in threshold values were found.')
    if max_rpos_diff > stringency:
        print('ALERT! Detected a position percentage difference of: ', max_rpos_diff*100, '%')
        # then indicate where the differences are
    else:
        print('Excellent! No differences in position values were found.')

    # Now check CSV files
