import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- Define parameters ---
subjects = ['S22','S29','S38','S40','S41','S42','S43','S46','S47',
             'S48','S49R','S50','S52','S53','S54','S55','S56','S57']

stdrel_vals = [0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
actr_vals = [50, 100, 150]

# --- Initialize storage for group data ---
mp_all = []
tp_all = []
pos_all = []

# --- Loop across subjects ---
for subj in subjects:
    mp_matrix = pd.DataFrame(np.nan, index=stdrel_vals, columns=actr_vals)
    tp_matrix = pd.DataFrame(np.nan, index=stdrel_vals, columns=actr_vals)
    pos_matrix = pd.DataFrame(np.nan, index=stdrel_vals, columns=actr_vals)

    for std in stdrel_vals:
        for act in actr_vals:
            folder_name = f"/Users/nicoletomassi/PycharmProjects/implant_model_v2_1/INV_OUTPUT/{subj}/combined_FIT_new_R250_std_{std}_act_{act}"
            file_path = os.path.join(folder_name, "summary_inverse_fit_results.csv")

            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path)

                    mp_error = df.loc[0, "Mean MP Threshold Error"]
                    tp_error = df.loc[0, "Mean TP Threshold Error"]
                    pos_error = df.loc[0, "Position Error"]

                    mp_matrix.loc[std, act] = mp_error
                    tp_matrix.loc[std, act] = tp_error
                    pos_matrix.loc[std, act] = pos_error

                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
            else:
                print(f"File not found: {file_path}")

    # Store subject matrices
    mp_all.append(mp_matrix)
    tp_all.append(tp_matrix)
    pos_all.append(pos_matrix)

# --- Compute group averages (element-wise) ---
mp_avg = sum(mp_all) / len(mp_all)
tp_avg = sum(tp_all) / len(tp_all)
pos_avg = sum(pos_all) / len(pos_all)

# --- Choose which metric to plot ---
metric_to_plot = tp_avg  # Change to tp_avg or pos_avg if desired

# --- Plot heatmap ---
plt.figure(figsize=(8, 6))
sns.set(font_scale=1.5)  # Scales font sizes globally (roughly ~18 pt)
ax = sns.heatmap(
    metric_to_plot,
    xticklabels=actr_vals,
    yticklabels=stdrel_vals,
    cmap='viridis',
    annot=True,
    fmt=".2f",
    cbar_kws={'label': 'Error'}
)

# --- Label formatting ---
plt.xlabel("actr", fontsize=18)
plt.ylabel("relstd", fontsize=18)
plt.title("Average MP Error Across Participants", fontsize=18)

ax.tick_params(axis='both', which='major', labelsize=16)
ax.figure.axes[-1].yaxis.label.set_size(18)  # Colorbar label

plt.tight_layout()
plt.show()
