import numpy as np
import pandas as pd

def read_subject_data(subject_name):
    file_path ='/Users/nicoletomassi/PycharmProjects/implant_model_v2_1/Subject_data.xlsx'
    sheet_name=subject_name
    columns=["Espace", "CT", "MP (dB)", "TP (dB)", "M", "Clinical_T", "M_levels uA","M_levels_dB"]

    reading = pd.read_excel(file_path, sheet_name=sheet_name, usecols=columns)
    data_matrix={col: reading[col].to_numpy() for col in reading.columns}

    return [data_matrix]