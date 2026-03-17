# script to cycle through a set of parameters to try to optimize fit results

import common_params as cp
import FwdModel4 as fwd
import FwdModel4_2D as fwd2D
import InverseModelCombined as inv_mod
import numpy as np
import csv
import os

# params to vary
# external resistivity, stdrel, number of active neurons, minimization approach
## planning resistivities of 70, 125, 250, 500, 1000
res_ext = 250.0  # DOUBLE CHECK this in COMMON PARAMS. should double-check that this matches the field table being used
# specify some fixed params

# start by cycling through stdrel, # of active neurons
# stdrel_vals = [0.125, 0.2, 0.5, 1, 2, 4, 8, 16, 32]
# thrtarg_vals = [1, 2, 5, 10, 25, 50, 100, 200, 500, 1000]
#stdrel_vals = [0.5, 1, 2, 4, 8]
stdrel_vals =  [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0] #[
thrtarg_vals = [25, 50, 100, 500, 1000]
actr_vals=[50, 100, 150]

# stdrel_vals = [2]
# thrtarg_vals = [100]


# Set up variables to store output
n_std = len(stdrel_vals)
n_thr = len(thrtarg_vals)
n_actr=len(actr_vals)
#mp_err_summ = np.zeros((n_std, n_thr))
#tp_err_summ = np.zeros((n_std, n_thr))

mp_err_summ = np.zeros((n_std, n_actr))
tp_err_summ = np.zeros((n_std, n_actr))

#dist_err_summ = np.zeros((n_std, n_thr))
#dist_corr_sum = np.zeros((n_std, n_thr))

dist_err_summ = np.zeros((n_std, n_actr))
dist_corr_sum = np.zeros((n_std, n_actr))

#dist_err_summ = np.zeros((n_std, n_thr))
#dist_corr_sum = np.zeros((n_std, n_thr))

diff_err_surv = np.zeros((n_std, n_actr))
dist_corr_surv = np.zeros((n_std, n_actr))

#n_corr_sig = np.zeros((n_std, n_thr))
n_corr_sig = np.zeros((n_std, n_actr))


for stdrel in stdrel_vals:  # loop on stdrel values
    espace = 1.1  # set espace to get 2D fwd model
    #for thrtarg in thrtarg_vals:  # loop on target #
    for act in actr_vals:
        new_dir_suffix = 'new_R%d' % res_ext + '_' + 'std_%.1f' % stdrel + '_act_%d' % act #'_thr_%d' % thrtarg

        FWD_OUT_PRFIX = 'FWD_OUTPUT/'
        FWDOUTPUTDIR = FWD_OUT_PRFIX + new_dir_suffix
        INV_OUT_PRFIX = 'INV_OUTPUT/'
        INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix
        # save file with key params
        param_file = 'surv_params.txt'  # param file to pass to the other model scripts
        tempdata = np.zeros(4)  # 4 values
        tempdata[0] = res_ext
        tempdata[1] = stdrel
        tempdata[2] = act#thrtarg
        tempdata[3] = espace
        with open(param_file, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',')
            for i, row in enumerate(tempdata):
                data_writer.writerow([tempdata[i]])
        data_file.close()

        # Run FwdModel
        fwd.fwd_model_4('survey')
        # Run FwdModel2D for espace = 1.1 mm
        espace = 1.1
        # save file with key params
        param_file = 'surv_params.txt'
        tempdata = np.zeros(4)  # 4 values
        tempdata[0] = res_ext
        tempdata[1] = stdrel
        tempdata[2] = act #thrtarg
        tempdata[3] = espace
        with open(param_file, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',')
            for i, row in enumerate(tempdata):
                data_writer.writerow([tempdata[i]])
        data_file.close()
        fwd2D.fwd_model_2d('survey')  # Note that
        # Prepare to run FwdModel2D for espace = 0.85 mm

        espace = 0.85
        # save file with key params
        param_file = 'surv_params.txt'
        tempdata = np.zeros(4)  # 3 values
        tempdata[0] = res_ext
        tempdata[1] = stdrel
        tempdata[2] = act #thrtarg
        tempdata[3] = espace
        with open(param_file, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',')
            for i, row in enumerate(tempdata):
                data_writer.writerow([tempdata[i]])
        data_file.close()
        fwd2D.fwd_model_2d('survey')


        espace = 1.1  # Set espace back to 1.1 in surv_params.txt
        # save file with key params
        param_file = 'surv_params.txt'
        tempdata = np.zeros(4)  # 3 values
        tempdata[0] = res_ext
        tempdata[1] = stdrel
        tempdata[2] = act #thrtarg
        tempdata[3] = espace
        with open(param_file, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',')
            for i, row in enumerate(tempdata):
                data_writer.writerow([tempdata[i]])
        data_file.close()

        # Run inverse model
        inv_mod.inverse_model_combined('survey')
        print('back from inverse model')

        # Collect and collate results stats ( do this in a separate script)

# text from trial
        # TST_OUT_DIR = 'TSTOUT/'
        # if not os.path.exists(TST_OUT_DIR):
        #     os.mkdir(TST_OUT_DIR)
        # dir_base = 'TEST_DIR_NAME/'
        # orig_dir = TST_OUT_DIR + dir_base
        # # Run inverse model
        # os.mkdir(orig_dir)
        #
        # # Rename fwd and inverse output directories
        # new_dir_suffix = '_R%d' % res_ext + '_' + 'std_%.1f' % stdrel + '_thr_%d' % thrtarg
        # new_fwd_dir = orig_dir[0:-1] + new_dir_suffix
        # os.rename(orig_dir, new_fwd_dir)
###

        # Rename fwd and inverse output directories
        new_dir_suffix = 'combined_FIT_new_R%d' % res_ext + '_' + 'std_%.1f' % stdrel + '_act_%d' % act #'_thr_%d' % thrtarg
        # offset = len(cp.FWD_OUT_PRFIX)
        FWD_OUT_PRFIX = 'FWD_OUTPUT/'
        FWDOUTPUTDIR = FWD_OUT_PRFIX + new_dir_suffix
        INV_OUT_PRFIX = 'INV_OUTPUT/'
        INVOUTPUTDIR = INV_OUT_PRFIX + new_dir_suffix

        new_fwd_dir = FWDOUTPUTDIR + new_dir_suffix
        os.rename(cp.FWDOUTPUTDIR, FWDOUTPUTDIR)
        # offset = len(cp.INV_OUT_PRFIX)
        new_inv_dir = INVOUTPUTDIR + new_dir_suffix
        os.rename(cp.INVOUTPUTDIR, INVOUTPUTDIR)

## clean up by removing param_file
os.remove(param_file)

# Run a separate script to calculate mean & std for each subject across parameters
# Want to find out the best parameter combination, and also how robust each subject fit is
# to changes in parameter sets
