#  Script to compare two voltage tables to look for differences
# This could be useful when changing stringency parameters in voltage_calc.py
# both tables have to have the same dimensions

import pickle
import numpy as np
import matplotlib.pyplot as plt

table1 = 'v_tables/16June2024_MedResolution_Rext250_Rint250_nonans.dat'
table2 = 'v_tables/18May2024_MedResolution_Rext70_nonans.dat'

with open(table1, 'rb') as combined_data:
    data = pickle.load(combined_data)
combined_data.close()

fp1 = data[0]
v_vals1 = data[1]
act_vals1 = data[2]

with open(table2, 'rb') as combined_data:
    data = pickle.load(combined_data)
combined_data.close()

fp2 = data[0]
v_vals2 = data[1]
act_vals2 = data[2]

ew = 2  # edge width
act_vals_diff = np.subtract(act_vals2, act_vals1)
act_vals_diff_non_edge = np.subtract(act_vals2[ew:-ew, ew:-ew], act_vals1[ew:-ew, ew:-ew])
act_vals_reldiff = np.divide(act_vals_diff, act_vals2)
act_vals_reldiff_non_edge = np.divide(act_vals_diff_non_edge, act_vals2[ew:-ew, ew:-ew])

print('max of act_vals_diff: ', np.nanmax(np.abs((act_vals_diff[:]))))
print('max of act_vals_reldiff: ', np.nanmax(np.abs((act_vals_reldiff[:]))))
print('max of act_vals_diff_non_edge: ', np.nanmax(np.abs(act_vals_diff_non_edge[:])))
print('max of act_vals_reldiff_non_edge: ', np.nanmax(np.abs((act_vals_reldiff_non_edge[:]))))
