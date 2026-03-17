# get_thresholds.py version of 25 October 2023
# Original code by Steven Bierer
# Translated to python by David J. Perkel
import pickle
import csv
import numpy as np
import thr_function as thr_f
import cylinder3d_makeprofile as c3dm
from common_params import *  # import common values across all models
import get_thresholds as gt
import surv_full
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import read_subject_data



fit_precision = 0.01

#subject = scenarios[0]


def compare_neuron_count(subject):
  #  if mode == 'main':
     #   pass  # proceed as normal
    plotfigs='False'


    plotfigs_1='True'
    if "fieldTable" not in locals():
        with open(FIELDTABLE, 'rb') as combined_data:
            data = pickle.load(combined_data)
            combined_data.close()
            # (output is 'fieldTable' and 'fieldParams')
            fp = data[0]
            # Temp fixup
            fp['zEval'] = np.array(fp['zEval'])
            act_vals = data[2]  # the data[1] has voltage values, which we are ignoring here
            simParams['grid']['table'] = act_vals
            GRID['table'] = act_vals
            COCHLEA['res1'] = fp['resInt'] * np.ones(NELEC)  # Note these valus do not match thos of Goldwyn et al., 2010
            COCHLEA['res2'] = fp['resExt'] * np.ones(NELEC)  # resistivities are in Ohm*cm (conversion to Ohm*mm occurs later)
            GRID['r'] = fp[
                'rspace']  # only 1 of the 3 cylindrical dimensions can be a vector (for CYLINDER3D_MAKEPROFILE)
            survVals = fit_survival#

            rposVals = fit_distance #CT_data #make sure to subtract from 1 if plotting
            MP_data=fit_MP
            TP_data=fit_TP

            hires = '_hi_res'
            nSurv = len(survVals)
            nRpos = len(rposVals)
            ##
            electrodes['rpos'] = np.zeros(NELEC)
            erb_empirical_array_M = []
            erb_empirical_array_T = []
            ELEC_MIDPOINT = GRID['z'][-1] / 2.0  # electrode array midpoint
            array_base = -(np.arange(NELEC - 1, -1, -1) * espace)
            array_mid = (array_base[0] + array_base[-1]) / 2.0
            electrodes['zpos'] = (ELEC_MIDPOINT - array_mid) + array_base
            electrodes['rpos']=rposVals
            ##
            # electrodes['zpos'] = ELEC_BASALPOS - (np.arange(NELEC - 1, -1, -1) * espace)
            #CHANNEL['timestamp'] = datetime.datetime.now()
            COCHLEA['radius'] = np.ones(NELEC) * fp['cylRadius']  # rpos is in mm, so be sure they fit inside the radius
            simParams['cochlea'] = COCHLEA
            simParams['electrodes'] = electrodes
            simParams['channel'] = CHANNEL
            simParams['grid'] = GRID
            avec = np.arange(0, 1.01, .01)  # create the neuron count to neuron spikes transformation
            rlvec = NEURONS['coef'] * (avec ** 2) + (1 - NEURONS['coef']) * avec
            rlvec = NEURONS['neur_per_clust'] * (rlvec ** NEURONS['power'])
            rlvltable = np.stack((avec, rlvec))  # start with the e-field(s) created above, but remove the current scaling
            n_sig = len(sigmaVals)
            #thr_sim_db = np.empty((nSurv, nRpos, n_sig))  # Array for threshold data for different elecs and sigma values
           # thr_sim_db[:] = np.nan
            #neuronact = np.empty((nSurv, nRpos, n_sig, 1, len(GRID['z'])))  # Only 1 electrode in this model
           # neuronact[:] = np.nan
            simParams['neurons']['rlvl'] = rlvltable
           # comparisons = np.zeros((nRpos, 15,2)) #change if upping the targets
            simParams['neurons']['nsurvival'] = \
                surv_full.surv_full(simParams['electrodes']['zpos'], survVals, simParams['grid']['z'])

            #    print('surv  ', j, ' of ', len(survVals))
               # #Get survival values for all 330 clusters from the 16 values at electrode positions
            #simParams['neurons']['nsurvival'] = \
                    #surv_full.surv_full(simParams['electrodes']['zpos'], survVals, simParams['grid']['z'])
            elec_vals = range(0, NELEC)
            n_neur_results=[]
            n_neur_offsets_results = []
            n_neur_offsets_M_db_results=[]
            n_neur_offsets_MP_results=[]#for j in elec_vals:  # Loop on stimulating electrodes TODO: get rid of this
            n_neur_results_500=[]
            n_neur_results_100=[]
            n_neur_results_500_TP = []
            posvals = np.arange(0, 33, 0.01)
            posvals -= posvals[-1] / 2.0
            fig, axs = plt.subplots(4, 4, figsize=(8, 6))
         #   fig2, axs2 = plt.subplots(4, 4, figsize=(8, 6))
            n_electrodes=16
            MP_vals=[]
            electrode_positions = (np.arange(n_electrodes) - (n_electrodes - 1) / 2) * espace
            for j, surv in enumerate(survVals):
                simParams['channel']['number'] = j
                simParams['channel']['sigma']=0
                aprofile = c3dm.cylinder3d_makeprofile(act_vals, fp, simParams)
                e_field = abs(aprofile)

                #for k in range(len(elec_vals)):


                    # Call function and get n_neur
                #n_neur, nvals, _aa = thr_f.thr_function(M_levels_uA[j], e_field, simParams)
               # MP_vals[j]=gt.get_thresholds(act_vals, fp, simParams)
                n_neur_offsets, nvals_offsets, _aa_offsets = thr_f.thr_function(M_levels_uA_offsets[j], e_field, simParams)
                n_neur_MP, nvals_MP, _aa_MP = thr_f.thr_function(MP_uA[j], e_field,simParams)
                area_M = np.trapz(nvals_offsets / 10, x=posvals)  # Area under the curve
                peak_M = np.max(nvals_offsets / 10)  # Peak response

                area_T = np.trapz(nvals_MP / 10, x=posvals)  # Area under the curve
                peak_T = np.max(nvals_MP / 10)  # Peak response
                erb_empirical_T = area_T / peak_T
                erb_empirical_array_T.append(erb_empirical_T)

                erb_empirical_M = area_M / peak_M
                erb_empirical_array_M.append(erb_empirical_M)

                normalized_nvals_M = (nvals_offsets / 10) / peak_M

                # Find indices where the signal crosses 0.5
                above_half_M = normalized_nvals_M >= 0.5
                crossings_M = np.where(np.diff(above_half_M.astype(int)) != 0)[0]

                if len(crossings_M) >= 2:
                    # Interpolate to get better precision on the crossing points
                    left_idx = crossings_M[0]
                    right_idx = crossings_M[-1]

                    # Linear interpolation for left crossing
                    x1, x2 = posvals[left_idx], posvals[left_idx + 1]
                    y1, y2 = normalized_nvals_M[left_idx], normalized_nvals_M[left_idx + 1]
                    left_cross = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

                    # Linear interpolation for right crossing
                    x1, x2 = posvals[right_idx], posvals[right_idx + 1]
                    y1, y2 = normalized_nvals_M[right_idx], normalized_nvals_M[right_idx + 1]
                    right_cross = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

                    half_peak_width_M = right_cross - left_cross
                else:
                    half_peak_width_M = np.nan  # Could not determine width

                normalized_nvals_T = (nvals_MP / 10) / peak_T

                # Find indices where the signal crosses 0.5
                above_half_T = normalized_nvals_T >= 0.5
                crossings_T = np.where(np.diff(above_half_T.astype(int)) != 0)[0]

                if len(crossings_T) >= 2:
                    # Interpolate to get better precision on the crossing points
                    left_idx = crossings_T[0]
                    right_idx = crossings_T[-1]

                    # Linear interpolation for left crossing
                    x1, x2 = posvals[left_idx], posvals[left_idx + 1]
                    y1, y2 = normalized_nvals_T[left_idx], normalized_nvals_T[left_idx + 1]
                    left_cross = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

                    # Linear interpolation for right crossing
                    x1, x2 = posvals[right_idx], posvals[right_idx + 1]
                    y1, y2 = normalized_nvals_T[right_idx], normalized_nvals_T[right_idx + 1]
                    right_cross = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)

                    half_peak_width_T = right_cross - left_cross
                else:
                    half_peak_width_T = np.nan  # Could not determine width



                #n_neur_M_dB, nvals_M_dB, _aa_M_dB = thr_f.thr_function(M_levels_dB_offset[j], e_field,simParams)
                #n_neur_500, nvals_500,_aa_500 = thr_f.thr_function(MP_data[j], e_field,simParams)
                #n_neur_500_TP, nvals_500_TP, _aa_500_TP=thr_f.thr_function(TP_data[j], e_field, simParams)
                #n_neur_100, nvals_100, _aa_100 = thr_f.thr_function(MP_data[j], e_field, simParams)
              #  n_neur_500_TP, nvals_500_TP, _aa_500_TP = thr_f.thr_function(TP_data[j], e_field, simParams)

                # n_neur_results[j]=n_neur

                    # Store results
                #n_neur_results.append((j, n_neur))
                n_neur_offsets_results.append((j, n_neur_offsets))
                n_neur_offsets_MP_results.append((j, n_neur_MP))
                #n_neur_offsets_M_db_results.append((j, n_neur_M_dB))
                #n_neur_results_500.append((j,n_neur_500))
                #n_neur_results_100.append((j, n_neur_500))
                #n_neur_results_500_TP.append((j, n_neur_500_TP))
                if plotfigs_1=='True':

                    if j < 4:
                        #label_text_1 = f"sim Thresholds (dB): : {MP_data[j]: .2f} \nSurvVals: {survVals[j]: .2f}" #\n {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        #label_text_1 = f"dist: {1-rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        label_text_1=f"dist: {1-rposVals[j]: .2f} \nSurv: {survVals[j]: .2f} \n#neur: {n_neur_offsets_results[j][1]: .2f}"
                    #    label_text_2 = f"T neur: {n_neur_offsets_MP_results[j][1]: .2f} \nERB: : {erb_empirical_T: .2f} \n WHP: {half_peak_width_T:.2f}"
                        axs[0,j].plot(posvals,nvals_offsets/10, label = label_text_1)
                        axs[0,j].plot(posvals, nvals_MP/10)
                        axs[0, j].axvline(electrode_positions[j], color='gray', linestyle='--', lw=1)

                        axs[0,j].set_title(f"Electrode {j + 1}")
                        #axs[0, j].set_ylim(0, .12)
                        axs[0,j].legend(fontsize=8)

                        #label_text_2 = f"dist: {1 - rposVals[j]:.2f} \nMP T (CU): {MP_clinical[j]: .2f} \nMP T (uA): {MP_dB[j]:.2f} \nMP T (dB): {MP_dB[j]:.2f} \nN count: {n_neur_offsets_MP_results[j][1]:.2f}"
                        #axs[0, j].plot(posvals, nvals_MP / 10,label = label_text_2)
                        #axs[0, j].set_title(f"Electrode {j + 1}")
                        #axs[0, j].set_ylim(0, .04)
                       # axs[0, j].legend(fontsize=6)

                    elif j > 3 and j < 8:
                        #label_text_1 = f"sim M dB: : {MP_data[j]: .2f} \nSurvVals: {survVals[j]: .2f}" #\n {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        # label_text_1 = f"dist: {1-rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        label_text_1 = f"fit dist: {1 - rposVals[j]: .2f} \nSurv: {survVals[j]: .2f} \n#neur: {n_neur_offsets_results[j][1]: .2f}"# \nERB: : {erb_empirical_M: .2f} \n WHP: {half_peak_width_M:.2f}"
                      #  label_text_2=f"T neur: {n_neur_offsets_MP_results[j][1]: .2f} \nERB: : {erb_empirical_T: .2f} \n WHP: {half_peak_width_T:.2f}"
                        axs[1, j-4].plot(posvals,nvals_offsets/10, label = label_text_1)
                        axs[1,j-4].plot(posvals, nvals_MP/10)
                        axs[1, j - 4].axvline(electrode_positions[j], color='gray', linestyle='--', lw=1)
                        #label_text = f"dist: {1 - rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][ 1]:.2f}"
                        #axs[1, j-4].plot(posvals, nvals_offsets/10, label = label_text)
                        axs[1,j-4].set_title(f"Electrode {j + 1}")
                        #axs[1, j - 4].set_ylim(0, .2)
                        axs[1, j-4].legend(fontsize=8)

                        #label_text_2 = f"dist: {1 - rposVals[j]:.2f} \nMP T (CU): {MP_clinical[j]: .2f} \nMP T (uA): {MP_dB[j]:.2f} \nMP T(dB): {MP_dB[j]:.2f} \nN count: {n_neur_offsets_MP_results[j][1]:.2f}"
                        #axs[1, j - 4].plot(posvals, nvals_MP / 10, label = label_text_2)
                        #axs[1, j - 4].set_title(f"Electrode {j + 1}")
                        #axs[1, j - 4].set_ylim(0, .04)
                        #axs[1, j - 4].legend(fontsize=6)

                    elif j > 7 and j<12:
                        #label_text_1 = f"sim M dB: : {MP_data[j]: .2f} \nSurvVals: {survVals[j]: .2f}" #label_text_1 = f"dist: {1 - rposVals[j]:.2f} \nsim T dB: : {MP_data[j]: .2f}"# \nSurvVals: {survVals[j]: .2f}"  # \n {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        # label_text_1 = f"dist: {1-rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        label_text_1 = f" fit dist: {1 - rposVals[j]: .2f} \nSurv: {survVals[j]: .2f} \n#neur: {n_neur_offsets_results[j][1]: .2f}"# \nERB: : {erb_empirical_M: .2f} \n WHP: {half_peak_width_M:.2f}"
                      #  label_text_2 = f"T neur: {n_neur_offsets_MP_results[j][1]: .2f} \nERB: : {erb_empirical_T: .2f} \n WHP: {half_peak_width_T:.2f}"
                        axs[2, j-8].plot(posvals,nvals_offsets/10, label = label_text_1)
                        axs[2,j-8].plot(posvals, nvals_MP/10)
                        axs[2, j - 8].axvline(electrode_positions[j], color='gray', linestyle='--', lw=1)
                        #label_text_1 = f"dist: {1 - rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        #axs[2, j-8].plot(posvals, nvals_offsets / 10, label=label_text_1)
                        axs[2,j - 8].set_title(f"Electrode {j + 1}")
                       # axs[2, j - 8].set_ylim(0, .12)
                        axs[2, j-8].legend(fontsize=8)

                        #label_text_2 = f"dist: {1 - rposVals[j]:.2f} \nMP T (CU): {MP_clinical[j]: .2f} \nMP T (uA): {MP_dB[j]:.2f} \nMP T (dB): {MP_dB[j]:.2f} \nN count: {n_neur_offsets_MP_results[j][1]:.2f}"
                        #axs[2, j - 8].plot(posvals, nvals_MP / 10, label = label_text_2)
                        #axs[2, j - 8].set_title(f"Electrode {j + 1}")
                        #axs2[2, j - 8].set_ylim(0, .04)
                        #axs[2, j - 8].legend(fontsize=6)

                    else:
                        label_text_1 = f"fit dist: {1 - rposVals[j]: .2f} \nSurv: {survVals[j]: .2f} \n#neur: {n_neur_offsets_results[j][1]: .2f}"# \nERB: : {erb_empirical_M: .2f} \n WHP: {half_peak_width_M:.2f}"
                      #  label_text_2 = f"T neur: {n_neur_offsets_MP_results[j][1]: .2f} \nERB: : {erb_empirical_T: .2f} \n WHP: {half_peak_width_T:.2f}"
                        #label_text_1 = f"sim M dB: : {MP_data[j]: .2f} \nSurvVals: {survVals[j]: .2f}" #label_text_1 = f"dist: {1 - rposVals[j]:.2f} \nsim T dB: {MP_data[j]: .2f}"# \nSurvVals: {survVals[j]: .2f}"  # \n {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        # label_text_1 = f"dist: {1-rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM(dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        axs[3, j-12].plot(posvals,nvals_offsets/10, label = label_text_1)
                        axs[3,j-12].plot(posvals, nvals_MP/10)
                        axs[3, j - 12].axvline(electrode_positions[j], color='gray', linestyle='--', lw=1)
                        #label_text_1 = f"dist: {1 - rposVals[j]:.2f} \nM(CU): {M_levels[j]: .2f} \nM (uA): {M_levels_uA_offsets[j]:.2f} \nM dB): {M_levels_dB_offset[j]:.2f} \nN count: {n_neur_offsets_results[j][1]:.2f}"
                        #axs[3, j - 12].plot(posvals, nvals_offsets / 10, label=label_text_1)
                        axs[3, j-12].set_title(f"Electrode {j + 1}")
                        #axs[3, j - 12].set_ylim(0,.12)
                        axs[3, j - 12].legend(fontsize=8)
                        #plt.tight_layout()

                        #label_text_2 = f"dist: {1 - rposVals[j]:.2f} \nMP T(CU): {MP_clinical[j]: .2f} \nMP T(uA): {MP_dB[j]:.2f} \nMP T(dB): {MP_dB[j]:.2f} \nN count: {n_neur_offsets_MP_results[j][1]:.2f}"
                        #axs[3, j - 12].plot(posvals, nvals_MP / 10, label = label_text_2)
                        #axs[3, j - 12].set_title(f"Electrode {j + 1}")
                       # axs[3, j - 12].set_ylim(0, .04)
                        #axs[3, j - 12].legend(fontsize=6)
                        plt.tight_layout()

                    fig.supxlabel('Longitudinal distance (mm)')
                    fig.supylabel('Fractional neuronal activation')
                    fig.tight_layout()


                    #plt.show()


        if plotfigs=='True':

            fig, axes = plt.subplots(1, 3, figsize=(10, 4))  # (rows, columns)
            for k in range(nRpos):
                if n_neur_offsets_results[k,1]==0:
                    continue
                else:
                #axes[0,0].plot(M_levels[k],  n_neur_results[k,1], marker ='o', linestyle='-',alpha =0.5, label=f'dist = {(1-rposVals[k]):.2f}')
                #axes[0, 0].axhline(y=100, color='r', linestyle='--', label='thresholds')
               # axes[0,1].plot(M_levels_dB[k], n_neur_results[k, 1], marker='o', linestyle='-', alpha=0.5)
                #axes[0, 1].axhline(y=100, color='r', linestyle='--')
                    axes[0].plot(M_levels_offsets[k], n_neur_offsets_results[k, 1], marker='o', linestyle='-', alpha=0.5)
                    axes[0].axhline(y=100, color='r', linestyle='--')
                    axes[2].plot(M_levels_dB_offset[k], n_neur_offsets_results[k, 1], marker='o', linestyle='-', alpha=0.5,label=f'dist = {(1-rposVals[k]):.2f}')
                    axes[2].axhline(y=100, color='r', linestyle='--')
                    axes[1].plot(M_levels_uA_offsets[k], n_neur_offsets_results[k, 1], marker='o', linestyle='-', alpha=0.5,
                                    label=f'dist = {(1 - rposVals[k]):.2f}')
                    axes[1].axhline(y=100, color='r', linestyle='--')
                    # Red dotted line at y = 100



       # axes[0,0].set_xlabel('M levels (CU)')
        #axes[0,0].set_ylabel('neuron count')
       # axes[0, 1].set_xlabel('M levels dB')
       # axes[0, 1].set_ylabel('neuron count')

        #axes[0].set_xlabel('M levels (CU)')
        #axes[0].set_ylabel('neuron count')
        #axes[2].set_xlabel('M levels (dB)')
        #axes[2].set_ylabel('neuron count')
        #axes[1].set_xlabel('M levels (uA)')
        #axes[1].set_ylabel('neuron count')
       # plt.grid(True)
      #  axes[0,0].set_title('S22')
        #axes[1, 0].set_title('S22-offsets')

        # Add a title centered above the top two plots
            fig.text(0.5, 0.97, subject , fontsize=12, ha='center')

        # Add a title centered above the bottom two plots
        #fig.text(0.5, 0.5, subject + "-offsets", fontsize=12, ha='center')

            plt.legend( title = 'dist',  bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
            plt.tight_layout()
            plt.show()



#    print(f"{n_neur}", f"{n_neur_offsets_results}")
    return [n_neur_offsets_MP_results, n_neur_offsets_results, survVals, rposVals] #probably also want to return dist and survival and M level db

#if __name__ == '__main__':
# compare_neuron_count('main')  # alternatives are 'main', 'gui' and "survey'

subjects = ['S22','S29','S38','S40', 'S41','S42','S43','S46','S47','S48','S49R','S50','S52','S53','S54','S55','S56','S57']#, 'S22''S29', 'S40', 'S43', 'S47', 'S49R', 'S50']
empty_list = []
full_neuron_results=np.zeros((16,1))
full_offset_neuron_results=np.zeros((16,1))
full_distances=np.zeros((16,1))
full_survival=np.zeros((16,1))


full_neuron_MPs=[]
M_neurons=[]


for s,subj in enumerate(subjects):
    subject_data = read_subject_data.read_subject_data(subj)

    datafile = "INV_OUTPUT/RE250_RI70std_0.6_thr_100/" + subj + "_fitResults_combined.csv"
    file = open(datafile)
    numlines = len(file.readlines())
    file.close()

    with open(datafile, newline='') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',')
        # Convert the data into a list of rows
        rows = [row for row in datareader]

    # Convert rows into a NumPy array (for easier column extraction)
    data_array = np.array(rows[1:17], dtype=float)  # Assuming data is numeric

    # Extract each column into a separate NumPy array
    ncol = data_array.shape[1]  # Get the number of columns

    # Create a list to hold each column as a separate NumPy array
    columns = [data_array[:, i] for i in range(ncol)]

    # If you want to assign each column to a separate variable:
    col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12 = columns

    offsets = col12
    fit_survival = col9
    fit_MP=col4
    fit_TP=col5
    fit_distance=col7
    CT_data = subject_data[0]['CT']
    MP_dB = subject_data[0]['MP (dB)']
    TP_dB = subject_data[0]['TP (dB)']
    M_levels = subject_data[0]['M']
    #M_levels_uA = subject_data[0]['M_levels uA']
    MP_clinical=subject_data[0]['Clinical_T']
    espace = subject_data[0]['Espace'][0]
    M_levels_dB = subject_data[0]['M_levels_dB']
    M_levels_dB_offset = M_levels_dB - offsets
    M_levels_uA_offsets = 10 ** (M_levels_dB_offset / 20)
    #M_levels_offsets = M_levels_uA_offsets * 1.2458
    ##account for offsets
    MP_uA=10**(fit_MP/20)

    neuron_results, offsets_neuron_results, survival, distance_fit = compare_neuron_count(subj)
    second_elements=[item[1] for item in neuron_results]
    neuron_array=np.array(second_elements)
    second_elements_1=[item[1] for item in offsets_neuron_results]
    M_neuron_array=np.array(second_elements_1)
    full_neuron_results=np.column_stack((full_neuron_results, neuron_array))
    full_offset_neuron_results=np.column_stack((full_offset_neuron_results, M_neuron_array))
    full_distances=np.column_stack((full_distances, distance_fit))
    full_survival=np.column_stack((full_survival, survival))

print(full_neuron_results)

#breakpoint
        #thr_sim_db[j, k, i] = tempvals[0].item()
                     #   nexttemp2 = tempvals[1]
                      #  if np.max(nexttemp2) == 0:
                       #     print("flat activation")
                       # neuronact[j, k, i, :, :] = nexttemp2
                        # neuroncount = np.sum(neuronact[j, k, i, :, :])  ## commented out because unused 24 May 2024


  #  targets= np.array([100,500, 1000])
   # for targ_count in targets:
      #  sim_params['neurons']['thrtarg'] = targ_count
     #   neuron_counts = gt.get_thresholds(act_vals,fp,sim_params) #TODO go through this
   # return
#
# def get_thresholds(field_table, field_params, sim_params):
#     # function that is called by both forward and inverse models
#     # Preallocate the arrays needed
#     nz = len(sim_params['grid']['z'])
#     if isinstance(sim_params['electrodes']['rpos'], list):
#         nelec = len(sim_params['electrodes']['rpos'])
#     elif not isinstance(sim_params['electrodes']['rpos'], np.ndarray):  # special case for 2D version of fwd model
#         nelec = 1
#     else:  # regular situation when called from forward or inverse model
#         nelec = len(sim_params['electrodes']['rpos'])
#     thresholds = np.empty(nelec)  # Array to hold threshold data for different stim electrodes and varied sigma values
#     thresholds[:] = np.nan
#     e_field = np.empty(nz)
#     e_field[:] = np.nan
#
#     if nelec == 1:
#         elec_vals = range(7, 8)  # use electrode 7 to be in middle of cochlear length
#         # set the position manually in the center of the array
#         espace = sim_params['electrodes']['zpos'][elec_vals[0]] - sim_params['electrodes']['zpos'][elec_vals[0] - 1]
#         neuron_midpoint = np.max(sim_params['grid']['z'])/2.0
#         # TODO -- can eliminate these lines?
#         # elec_start = sim_params['electrodes']['zpos'][0]
#         # elec_end = sim_params['electrodes']['zpos'][-1]
#         elec_pos = neuron_midpoint  # fixed position in center of neuron array
#         # elec_pos = (elec_start + elec_end)/2.0
#         sim_params['electrodes']['zpos'][elec_vals[0]] = elec_pos
#         sim_params['electrodes']['zpos'][elec_vals[0]-1] = elec_pos - espace
#         sim_params['electrodes']['zpos'][elec_vals[0]+1] = elec_pos + espace
#     else:
#         if sim_params['channel']['sigma'] == 0.0:  # monopolar
#             elec_vals = range(0, nelec)
#         else:
#             elec_vals = range(1, nelec - 1)  # Can't do tripolar stimulation at end electrodes
#
#     # print('electrode positions: ', sim_params['electrodes'])
#     nvalarray = []
#     nvals = 0
#     counter = -1
#     for j in elec_vals:  # Loop on stimulating electrodes
#         targets = np.array([100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000])
#         threshold_comparison = np.zeros((np.size(targets), 2))
#         for targ_count in targets:
#             counter = counter + 1
#             sim_params['neurons']['thrtarg'] = targ_count
#
#             sim_params['channel']['number'] = j
#             aprofile = c3dm.cylinder3d_makeprofile(field_table, field_params, sim_params)
#             # This is the biophysically correct behavior, to include sidelobes
#             e_field = abs(aprofile)  # uV^2/mm^2
#
#             # We can remove sidelobes by setting negative values to zero
#             # for kk in range(len(aprofile):
#             #     if aprofile[kk] < 0.0:
#             #         efield[kk] = 0.0
#             #     else:
#             #         efield[kk] = aprofile[kk]
#             # e_field = aprofile
#             # e_field[np.where(aprofile < 0.0)] = 0.0
#
#             # for kk in range(len(e_field)):
#             #     if aprofile[kk] < 0.0:
#             #         e_field[kk] = 0.0
#             #     else:
#             #         e_field[kk] = aprofile[kk]
#
#             # Use our own solving algorithm; this is clumsy but I couldn't find a solver in scipy.optimize that took
#             # advantage of the known monotonicity of the function
#             error = 20  # Ensure that the system starts with a large "error"
#             target = sim_params['neurons']['thrtarg'] #TODO make this go through a number of targets
#             nextstim = sim_params['channel']['current']
#             lastpos = nextstim
#             lastneg = 0.0
#             while np.abs(error) > fit_precision:
#                 [n_neur, nvals, _aa] = thr_f.thr_function(nextstim, e_field, sim_params)
#                 # returns nCount, NeuralProfile, ActiveProfile
#                 error = n_neur - target
#                 if error < -fit_precision:  # want to increase stim
#                     if nextstim == sim_params['channel']['current']:
#                         break
#                     lastneg = nextstim
#                     nextstim = nextstim + ((lastpos - nextstim)/2.0)
#                 elif error > fit_precision:  # decrease stimulus
#                     lastpos = nextstim
#                     nextstim = nextstim - ((nextstim - lastneg)/2.0)
#             if nelec == 1:
#                 thresholds[0] = nextstim
#             else:
#                 thresholds[j] = nextstim

            # Here call thrFunction to get the neural profile at threshold
            # if len(thresholds) > 1:
            #     profileval = thr_f.thr_function(thresholds[j], e_field, sim_params)[1]
            #     nvalarray.append(profileval)
            #
            # thresholds = 20 * np.log10(thresholds)  # return thresholds in dB
            # ncount = np.sum(np.asarray(nvals))
     #       threshold_comparison[counter, 0]=thresholds
    #        threshold_comparison[counter, 1]= ncount

   # return threshold_comparison


