# load_fwd_csv_data.py
#  Utility function to load threshold, CT, and potentially other data from a
#  spreadsheet.
#  Translated to Python 3 by David Perkel December 2020

import csv
import numpy as np
from common_params import *


def load_fwd_csv_data(loadfile, mode):

    thr_data = {'thrmp_db': np.zeros(NELEC), 'thrmp': np.zeros(NELEC), 'thrtp_db': np.zeros(NELEC),
                'thrtp': np.zeros(NELEC), 'thrtp_sigma': 0.9, 'mp_mcl_db': np.zeros(NELEC)}
    ct_data = {'stdiameter': [], 'scala': [], 'elecdist': [], 'espace': 1.1, 'type': [], 'insrt_base': [],
               'insert_apex': []}

    radius = 1.0

    # Load the data
    if mode == 'mp_tp':
        # file format
        thr_mp_col = 3
        thr_tp_col = 4
        thisrow = 0
        with open(loadfile, mode='r') as data_file:
            data_reader = csv.reader(data_file, delimiter=',', quotechar='"')
            for row in data_reader:
                if row[thr_mp_col] == 'nan':
                    thr_data['thrmp_db'][thisrow] = np.nan
                else:
                    aaa = row[thr_mp_col]
                    bbb = np.fromstring(aaa, sep=' ')
                    # ccc = float(bbb)
                    thr_data['thrmp_db'][thisrow] = float(bbb[0])

                thr_data['thrmp'][thisrow] = np.power(10.0, (thr_data['thrmp_db'][-1] / 20))

                if row[thr_tp_col] == 'nan':
                    thr_data['thrtp_db'][thisrow] = np.nan
                else:
                    aaa = row[thr_tp_col]
                    bbb = np.fromstring(aaa, sep=' ')
                    thr_data['thrtp_db'][thisrow] = float(bbb[0])

                thr_data['thrtp'][thisrow] = np.power(10.0, (thr_data['thrtp_db'][-1] / 20))

            # thr_data['thrmp_db'] = np.array(thr_data['thrmp_db'])
            # thr_data['thrtp_db'] = np.array(thr_data['thrtp_db'])
                thisrow += 1

    elif mode == 'mp_mcl':
        mp_mcl_col = 3
        with open(loadfile, mode='r') as data_file:
            data_reader = csv.reader(data_file, delimiter=',', quotechar='"')
            for row in data_reader:
                if row[mp_mcl_col] == 'nan':
                    thr_data['mp_mcl_db'].append(np.nan)
                else:
                    thr_data['mp_mcl_db'].append(float(np.fromstring(row[mp_mcl_col], sep=' ')))
    elif mode == 'ecap':
        pass

    n_elec = len(thr_data['thrmp_db'])
    ct_data['stdiameter'] = radius * 2.0 * (np.zeros(n_elec) + 1.0)
    return [thr_data, ct_data]
