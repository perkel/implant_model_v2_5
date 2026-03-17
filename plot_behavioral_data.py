import matplotlib.pyplot as plt
import pandas as pd

subjects = ['S22','S29','S38','S40', 'S41','S42','S43','S46','S47','S48','S49R','S50','S52','S53','S54','S55','S56','S57']
file_name="behavioral_data.xlsx"

n_subjects=len(subjects)
# Load each sheet into a DataFrame, then convert to NumPy arrays
M_data = pd.read_excel(file_name, sheet_name="M_level", header=None).to_numpy()
MP_data = pd.read_excel(file_name, sheet_name="MP", header=None).to_numpy()
TP_data = pd.read_excel(file_name, sheet_name="TP", header=None).to_numpy()


fig, axes = plt.subplots(3, 6, figsize=(18, 9), sharex=True, sharey=True)
axes = axes.flatten()  # flatten into 1D for easy indexing

for subj in range(n_subjects):
    ax=axes[subj]
    ax.plot(M_data[:, subj], label="M_level", color="green")
    ax.plot(MP_data[:, subj], label="MP", color="blue")
    ax.plot(TP_data[:, subj], label="TP", color="red")
    ax.set_title(subjects[subj], fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.4)
# Add global title
fig.suptitle("Behavioral Data Across Subjects", fontsize=16, weight="bold", y=1.02)

# Add global x and y labels
fig.supxlabel("Electrodes", fontsize=14)
fig.supylabel("dB", fontsize=14)

# Add a single legend outside the subplots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper left", ncol=3, fontsize=12)

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.suptitle("Behavioral Data Across Subjects", fontsize=16, weight="bold")

plt.show()