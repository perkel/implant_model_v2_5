# Sanity check to compare 2D maps to see how much they change for different parameters or espace values

import matplotlib.pyplot as plt
from common_params import *
import numpy as np
import csv

f1 = '2D_0_5_STDR0_6_espace_110.csv'
f2 = '2D_0_9_STDR0_6_espace_110.csv'

# Load monopolar data
datafile = FWDOUTPUTDIR + f1
file = open(datafile)
numlines = len(file.readlines())
file.close()

with open(datafile, newline='') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',')
    ncol = len(next(datareader))
    csvfile.seek(0)
    f1_thr = np.empty([numlines, ncol])
    for i, row in enumerate(datareader):
        # Do the parsing
        f1_thr[i, :] = row

datafile = FWDOUTPUTDIR + f2
file = open(datafile)
numlines = len(file.readlines())
file.close()

with open(datafile, newline='') as csvfile:
    datareader = csv.reader(csvfile, delimiter=',')
    ncol = len(next(datareader))
    csvfile.seek(0)
    f2_thr = np.empty([numlines, ncol])
    for i, row in enumerate(datareader):
        # Do the parsing
        f2_thr[i, :] = row

diff_thr = np.abs(np.subtract(f2_thr, f1_thr))
frac_diff_thr = np.divide(diff_thr, f1_thr)

print('max frac diff: ', np.max(frac_diff_thr))

where = np.argmax(frac_diff_thr, keepdims = True)
print('max is at: ', where)

