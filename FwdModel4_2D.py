# Forward model 2D version

# This 2D version uses a single electrode and scans through an arbitrary set of values for survival values and
# electrode position
# It is critical to note that the parameter values do not to be tied to the # of electrodes in an implant
# David J. Perkel 7 March 2019
# Translated to python 8 March 2020
# Reorganizing the code to be more modular and to allow sharing of code
# with inverse model, for simplicity and reproducibility
# This version runs through a set of 'scenarios', each with a specific set
# of electrode radial positions and neuron survival values. For each, it
# calculates predicted threshold for monopolar and partial tripolar
# stimulation, typically sigma = 0 and 0.9.
# Re-written 12 August 2021 to simplify data structures into dicts

import surv_full
import get_thresholds as gt
import datetime
import os
import csv
import matplotlib.pyplot as plt
import pickle
import calc_ecap as c_ecap
from common_params import *
import figs4_5_2D_contour
import fig3_neuron_activation


def fwd_model_2d(mode):

    n_stim = 4000  # number of stim levels for ecap calculation
    # TODO:  should loop on espace values. Looping could then be removed from the param_survey script
    if mode == 'main':
        espace_range = [0.85, 1.1]

        pass  # proceed as normal
    elif mode == 'survey':
        # load survey params
        param_file = 'surv_params.txt'
        tempdata = np.zeros(4)  # 4 values
        if os.path.exists(param_file):
            with open(param_file, newline='') as csvfile:
                datareader = csv.reader(csvfile, delimiter=',')
                csvfile.seek(0)
                for i, row in enumerate(datareader):
                    # Do the parsing
                    tempdata[i] = row[0]

        #  res_ext = tempdata[0]  # this variable is not used in this context.
        NEURONS['act_stdrel'] = tempdata[1]
        NEURONS['thrtarg'] = 100
        NEURONS['act_ctr'] = tempdata[2]
        espace_range = tempdata[3]
    else:  # should not happen
        print('fwd_model called with unrecognized mode: ', mode)
        exit()

    # We depend on a voltage and activation tables calculated using voltage_calc.py and saved as a .dat file
    with open(FIELDTABLE, 'rb') as combined_data:
        data = pickle.load(combined_data)
        combined_data.close()

    # Unpack data from the field data file
    fp = data[0]
    #  vVals = data[1]  # Voltage values are loaded from the file but are not needed
    act_vals = data[2]
    # Convert zEval to a numpy array
    fp['zEval'] = np.array(fp['zEval'])

    GRID['table'] = act_vals

    COCHLEA['res1'] = fp['resInt'] * np.ones(NELEC)  # Note these valus do not match thos of Goldwyn et al., 2010
    COCHLEA['res2'] = fp['resExt'] * np.ones(NELEC)  # resistivities are in Ohm*cm (conversion to Ohm*mm occurs later)
    GRID['r'] = fp['rspace']  # only 1 of the 3 cylindrical dimensions can be a vector (for CYLINDER3D_MAKEPROFILE)

    # The lines below are not used, but represent infrastructure for nonlinear loudness growth curves or
    # ecap amplitude growth
    # coef = 0  # convex: -1 | 0.4, 0.9 ; linear: 0 | 1; concave: +1 | 1.0, 1.8
    # POW = 1
    # SIDELOBE = 1
    ################################

    if_plot = False  # Whether to plot the results
    survVals = np.arange(0.04, 0.97, 0.02)
    rposVals = np.arange(-0.95, 0.96, 0.02)
    # survVals = np.arange(0.04, 0.97, 0.2)  # fewer entries for faster debugging
    # rposVals = np.arange(-0.95, 0.96, 0.2)
    hires = '_hi_res'
    nSurv = len(survVals)
    nRpos = len(rposVals)

    # set up filename
    descrip = 'surv_%.2f' % np.min(survVals) + '_%.2f' % np.max(survVals) + "_rpos_%.2f" %\
              np.min(rposVals) + '_%.2f' % np.max(rposVals) + hires

    if not os.path.isdir(FWDOUTPUTDIR):
        os.mkdir(FWDOUTPUTDIR)


    for this_es in espace_range:
        if this_es == 0.85:
            e_txt = '085'
        elif this_es == 1.1:
            e_txt = '110'
        else:
            e_txt = 'xxx'
        es_text = '_espace_' + e_txt

        outfile = FWDOUTPUTDIR + 'FwdModelOutput_' + descrip + es_text + '.csv'

        # Additional setup
        COCHLEA['timestamp'] = datetime.datetime.now()
        electrodes['timestamp'] = datetime.datetime.now()
        electrodes['rpos'] = np.zeros(NELEC)
        ELEC_MIDPOINT = GRID['z'][-1] / 2.0  # electrode array midpoint
        array_base = -(np.arange(NELEC - 1, -1, -1) * this_es)
        array_mid = (array_base[0] + array_base[-1]) / 2.0
        electrodes['zpos'] = (ELEC_MIDPOINT - array_mid) + array_base

        CHANNEL['timestamp'] = datetime.datetime.now()
        COCHLEA['radius'] = np.ones(NELEC) * fp['cylRadius']  # rpos is in mm, so be sure they fit inside the radius

        # Construct the simParams structure
        simParams['cochlea'] = COCHLEA
        simParams['electrodes'] = electrodes
        simParams['channel'] = CHANNEL
        simParams['grid'] = GRID

        # Example of neural activation at threshold (variants of Goldwyn paper) ; using choices for channel, etc.,
        # made above. Keep in mind that the activation sensitivity will be FIXED, as will the number of neurons required
        # for threshold. Therefore, the final current to achieve theshold will vary according to the simple minimization
        # routine. Also note this works much faster when the field calculations are performed with a look-up table.
        avec = np.arange(0, 1.01, .01)  # create the neuron count to neuron spikes transformation
        rlvec = NEURONS['coef'] * (avec ** 2) + (1 - NEURONS['coef']) * avec
        rlvec = NEURONS['neur_per_clust'] * (rlvec ** NEURONS['power'])
        rlvltable = np.stack((avec, rlvec))  # start with the e-field(s) created above, but remove the current scaling

        n_sig = len(sigmaVals)
        # n_behavior = len(behaviorVals)
        #  how many entries in the neuronact 3rd dimension? One each for mp & tp threshold plus 1 for mp_mcl
        ndim = 0
        for entry in all_measures:
            if entry == 'mp':
                ndim += 1
            if entry == 'tp':
                ndim += 1
            if entry == 'mp_mcl':
                ndim += 1
        thr_sim_db = np.empty((nSurv, nRpos, n_sig))  # Array for threshold data for different elecs and sigma values
        thr_sim_db[:] = np.nan
        neuronact = np.empty((nSurv, nRpos, ndim, 1, len(GRID['z'])))  # Only 1 electrode in this model
        neuronact[:] = np.nan
        # neuroncount=[]

        # Get survival values for all neuron clusters from the 16 values at electrode positions.
        simParams['neurons']['rlvl'] = rlvltable
        m_level_db = np.empty((nSurv, nRpos))  # Array for threshold data for different stim elecs and sigma values
        m_level_db[:] = np.nan
        ecap_all = np.zeros((nSurv, nRpos, n_stim, 2))

        # Sanity check. Could add other sanity checks here
        # if any(simParams.grid.r < 1):
        # raise('Ending script. One or more evaluation points are inside cylinder; not appropriate for neural activation.')

        # Determine threshold for each value of sigma, looping through survival values and radial positions
        for j, surv in enumerate(survVals):
            print('surv  ', j, ' of ', len(survVals))
            # Get survival values for all 330 clusters from the 16 values at electrode positions
            simParams['neurons']['nsurvival'] =\
                surv_full.surv_full(simParams['electrodes']['zpos'], np.ones(NELEC)*surv, simParams['grid']['z'])

            for k, rpos in enumerate(rposVals):
                simParams['electrodes']['rpos'] = rpos

                if 'mp' or 'tp' in all_measures:
                    # Calculate momopolar and tripolar thresholds
                    for i in range(0, n_sig):  # number of sigma values to test
                        simParams['channel']['sigma'] = sigmaVals[i]
                        tempvals = gt.get_thresholds(act_vals, fp, simParams)
                        thr_sim_db[j, k, i] = tempvals[0].item()
                        nexttemp2 = tempvals[1]
                        if np.max(nexttemp2) == 0:
                            print("flat activation")
                        neuronact[j, k, i, :, :] = nexttemp2

                if 'mp_mcl' in all_measures:

                    simParams['channel']['sigma'] = 0  # sigmaVals[i]
                    simParams['channel']['behavior'] = 1  # search for m_level, not threshold
                    tempvals = gt.get_thresholds(act_vals, fp, simParams)
                    m_level_db[j, k] = tempvals[0].item()
                    nexttemp2 = tempvals[1]
                    if np.max(nexttemp2) == 0:
                        print("flat activation")
                    neuronact[j, k, i, :, :] = nexttemp2
                     # neuroncount = np.sum(neuronact[j, k, i, :, :])  ## commented out because unused 24 May 2024

                if 'ecap' in all_measures:
                    # calculate ecap growth curve
                    ec_thresh = 300  # neurons
                    # Calculate ecap growth curve for stim intensities within the range just defined
                    ecap_temp = c_ecap.calc_ecap(act_vals, fp, simParams, ec_thresh, thr_sim_db[:, i], n_stim)
                    ecap_all[j, k, :, :] = ecap_temp


        # Write results to a CSV file
        header1 = ['Monopolar thresholds', 'rpos values in columns']
        header2 = ['Survival ']
        for i, rpos in enumerate(rposVals):
            header2 += str(rpos) + ', '

        with open(outfile, mode='w') as data_file:
            data_writer = csv.writer(data_file, delimiter=',', quotechar='"')
            if 'mp' in all_measures:
                data_writer.writerow(header1)
                data_writer.writerow(header2)
                for row, s_val in enumerate(survVals):
                    data_writer.writerow([row, s_val, thr_sim_db[row, 0]])

            if 'tp' in all_measures:
                header1 = 'Tripolar thresholds: rpos values in columns; sigma = ' + str(sigmaVals[-1])
                data_writer.writerow(header1)
                data_writer.writerow(header2)
                for row, s_val in enumerate(survVals):
                    data_writer.writerow([row, s_val, thr_sim_db[row, 1]])

            if 'mp_mcl' in all_measures:
                header1 = 'Most comfortable levels: rpos values in columns; sigma = 0.0'
                data_writer.writerow(header1)
                data_writer.writerow(header2)
                for row, s_val in enumerate(survVals):
                    data_writer.writerow([row, s_val, m_level_db[row]])

        data_file.close()

        if 'mp' in all_measures:
            np.savetxt(FWDOUTPUTDIR + 'Monopolar_2D_' + STD_TEXT + es_text + '.csv', thr_sim_db[:, :, 0], delimiter=',')
        if 'tp' in all_measures:
            np.savetxt(FWDOUTPUTDIR + 'Tripolar_09_2D_' + STD_TEXT + es_text + '.csv', thr_sim_db[:, :, 1], delimiter=',')
        if 'mp_mcl' in all_measures:
            np.savetxt(FWDOUTPUTDIR + 'MP_mcl_2D_' + STD_TEXT + es_text + '.csv', m_level_db[:, :], delimiter=',')

        spname = FWDOUTPUTDIR + 'simParams' + descrip + es_text + '.pickle'
        with open(spname, 'wb') as f:
            pickle.dump(simParams, f, pickle.HIGHEST_PROTOCOL)
            f.close()
        print('saved: ', spname)
        # Note that this is saving only the last simParams structure from the loops on sigma and in get_thresholds.

        # display min and max threshold values
        print('min MP thr:  ', np.min(thr_sim_db[:, :, 0]))
        print('max MP thr:  ', np.max(thr_sim_db[:, :, 0]))
        print('min TP thr:  ', np.min(thr_sim_db[:, :, 1]))
        print('max TP thr:  ', np.max(thr_sim_db[:, :, 1]))

        # Save neuron activation data into a binary file
        np.savez(FWDOUTPUTDIR + 'neuronact_' + STD_TEXT + es_text, survVals, rposVals, neuronact)

        #  Save ecap data to a binary file
        np.savez(FWDOUTPUTDIR + 'ecap_' + STD_TEXT + es_text, survVals, rposVals, ecap_all)

    # Plot the results
    if if_plot:
        figs4_5_2D_contour.fig_2D_contour()
        fig3_neuron_activation.fig3_neuron_activation()
        plt.show()


if __name__ == '__main__':
    fwd_model_2d('main')  # mode is 'main' or 'survey'
