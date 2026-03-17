#script to generate Figure 2 of Modeling M levels -Nicole Tomassi 12/1/2025

import pickle
import datetime
import os
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#  These are local to this project
import set_scenario as s_scen
import surv_full
import get_thresholds as gt
from common_params import *
import thr_function as th
import cylinder3d_makeprofile as c3dm

plt.rcParams.update({
    'font.size': 18,          # base font size
    'axes.titlesize': 18,     # title font
    'axes.labelsize': 18,     # x/y label font
    'xtick.labelsize': 16,    # tick labels
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

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
    NEURONS['act_stdrel'] = tempdata[1]
    NEURONS['thrtarg'] = tempdata[2]
    espace = tempdata[3]  # should be overridden by scenario/subject

# We depend on voltage and activation tables calculated using
# voltage_calc.py and saved as a .dat file. the file is specified in common_params.py
with open(FIELDTABLE, 'rb') as combined_data:
    data = pickle.load(combined_data)
    combined_data.close()

fp = data[0]
fp['zEval'] = np.array(fp['zEval'])
#  v_vals = data[1] #  voltage values not used
act_vals = data[2]
GRID['table'] = act_vals

COCHLEA['res1'] = fp['resInt'] * np.ones(NELEC)  # Note these values do not match those of Goldwyn et al., 2010
COCHLEA['res2'] = fp['resExt'] * np.ones(NELEC)  # resistivities are in Ohm*cm (conversion to Ohm*mm occurs later)
GRID['r'] = fp['rspace']  # only 1 of the 3 cylindrical dimensions can be a vector (for CYLINDER3D_MAKEPROFILE)

if_plot = True  # Whether to plot the results
n_sig = len(sigmaVals)


for scenario in scenarios:
    # if this scenario is a subject, set use_forward_model to be false
    first_let = scenario[0]
    if (first_let == 'A' or first_let == 'S') and scenario[1:3].isnumeric():
        print('This scenario, ', scenario, ' appears to be for a subject, not a forward model scenario. Skipping.')
        #surv_vals=

    [surv_vals, electrodes['rpos'], espace] = s_scen.set_scenario(scenario, NELEC)
    #electrodes['rpos']=electrodes['rpos']*10
    electrodes['rpos']=electrodes['rpos']*10
    elec_midpoint = GRID['z'][-1]/2.0  # electrode array midpoint
    array_base = -(np.arange(NELEC - 1, -1, -1) * espace)
    array_mid = (array_base[0] + array_base[-1])/2.0
    # electrodes['zpos'] = ELEC_BASALPOS - (np.arange(NELEC - 1, -1, -1) * espace)

    electrodes['zpos'] = (elec_midpoint - array_mid) + array_base

    if not os.path.isdir(FWDOUTPUTDIR):
        os.makedirs(FWDOUTPUTDIR)

    outfile = FWDOUTPUTDIR + 'FwdModelOutput_' + scenario + '.csv'


    # Additional setup
    RUN_INFO['scenario'] = scenario
    RUN_INFO['run_time'] = datetime.datetime.now()
    COCHLEA['radius'] = np.ones(NELEC) * fp['cylRadius']  # note that '.rpos' is in mm, must fit inside the radius

    # Construct the simParams structure
    simParams['cochlea'] = COCHLEA
    simParams['electrodes'] = electrodes
    simParams['channel'] = CHANNEL
    simParams['grid'] = GRID
    simParams['run_info'] = RUN_INFO
    simParams['neurons']['thrtarg'] = 100
    simParams['channel']['number'] = 7

    # Example of neural activation at threshold (variants of Goldwyn paper) ; using choices for channel, etc.
    # made above.
    # Keep in mind that the activation sensitivity will be FIXED, as will the number of neurons required
    # for threshold. Therefore, the final current to achieve threshold will vary according to the simple
    # minimization routine.
    avec = np.arange(0, 1.01, .01)  # create the neuron count to neuron spikes transformation
    rlvec = NEURONS['coef'] * (avec ** 2) + (1 - NEURONS['coef']) * avec
    rlvec = NEURONS['neur_per_clust'] * (rlvec ** NEURONS['power'])
    rlvltable = np.stack((avec, rlvec))  # start with the e-field(s) created above, but remove the current scaling

    # Specify which variables to vary and set up those arrays
    thr_sim_db = np.empty((NELEC, n_sig))  # Array for threshold data for different stim elecs and diff sigma values
    thr_sim_db[:] = np.nan

    # Get survival values for all 330 clusters from the 16 values at electrode positions.
    NEURONS['nsurvival'] = surv_full.surv_full(simParams['electrodes']['zpos'], surv_vals, simParams['grid']['z'])
    NEURONS['rlvl'] = rlvltable
    simParams['neurons'] = NEURONS
    #  Sanity check. Could add other sanity checks here
    if simParams['grid']['r'] < 1:
        raise ('Ending script. One or more evaluation points are inside cylinder;\
         not appropriate for neural activation.')

    # Determine threshold for each value of sigma

    MP_thresholds=np.linspace(20,90,50)
    simParams['channel']['sigma'] = 0

    testing_act_vals=[.4, 2]

    fig, axs = plt.subplots(1, 2, figsize=(8, 6))
    posvals = np.arange(0, 33, 0.01)
    posvals -= posvals[-1] / 2.0

    fig_rows = 0
    #plot the stimulus curve

    #panel 1
    CHANNEL['number']=7
    simParams['channel']=CHANNEL # j
    for act in testing_act_vals:
        NEURONS['act_stdrel']=act
        simParams['neurons'] = NEURONS
        nvalarray = []
        nvals_count = []

        elec_vals = range(7, 8)
        p = 0
        max_fraction = np.zeros(np.size(MP_thresholds))
        for j in MP_thresholds:
            CHANNEL['number'] = 7
            simParams['channel']=CHANNEL
            aprofile = c3dm.cylinder3d_makeprofile(act_vals, fp, simParams)
            # This is the biophysically correct behavior, to include sidelobes
            e_field = abs(aprofile)
            simParams['channel']['sigma'] = 0  # 0.9
            current = j
            uA_current = 10 ** (current / 20)
            n_neur, nvals, _aa = th.thr_function(uA_current, e_field, simParams)
            nvalarray.append(nvals)
            nvals_count.append(n_neur)

            activation = nvalarray[p]
            # Total number of activated neurons
            activation = activation / 10
            N_total = np.sum(activation)
            # Peak fractional activation and index
            max_fraction[p] = np.max(activation)
            p = p + 1
        axs[0].plot(MP_thresholds,max_fraction, color= 'black')
        axs[0].set_xlabel("dB")
        axs[0].set_title("Electrodes Far from Neurons")


        #plot where on the curve we are to stimulate 100, 500, and 1000
        neuron_targs=[100,500,1000]
        colors=['red','blue', 'green']
        color_idx =0
        for targ in neuron_targs:
            CHANNEL['number'] = 7
            simParams['channel']=CHANNEL
            simParams['channel']['number'] = 7  # j
            simParams['neurons']['thrtarg'] = targ
            aprofile = c3dm.cylinder3d_makeprofile(act_vals, fp, simParams)
            e_field = abs(aprofile)
            [thr_sim_db, neuron_vals] = gt.get_thresholds(act_vals, fp, simParams)

            # This is the biophysically correct behavior, to include sidelobes


            MP_thresholds_new = thr_sim_db
            f = interp1d(MP_thresholds, max_fraction, kind='linear')  # or 'cubic', 'quadratic', etc.
            # Interpolate at a new x value
            x_new = MP_thresholds_new[1]
            y_new = f(x_new)
            axs[0].plot(x_new, y_new, 'o', color=colors[color_idx])

            uA_current = 10 ** (x_new / 20)
            n_neur, nvals, _aa = th.thr_function(uA_current, e_field, simParams)
            nvalarray.append(nvals)
            axs[1].plot(posvals,nvals/10, color=colors[color_idx])
            #axs[1].set_ylim([0,0.3])
            axs[1].set_xlim([-5, 5])

            color_idx=color_idx+1
            if color_idx==3:
                axs[1].set_xlabel('dB')
                axs[1].set_xlabel('Longitudinal distance (mm)')
                axs[0].set_ylabel('Fractional Neuronal Activation')

        fig_rows=fig_rows+1

    plt.tight_layout()
    plt.show()
breakpoint()


      #  f = interp1d(MP_thresholds, max_fraction, kind='linear')  # or 'cubic', 'quadratic', etc.
        # Interpolate at a new x value
     #   x_new = fit_MP[k]
     #   y_new = f(x_new)
    #    plt.plot(x_new, y_new, 'o', color='red')

