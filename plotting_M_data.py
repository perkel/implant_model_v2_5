import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# Load your CSVs
M_level = pd.read_csv("M_levels_db.csv")     # shape (16, 18)
MP_data = pd.read_csv("MP_db.csv")     # shape (16, 18)
M_neurons=pd.read_csv("full_M_neuroncount.csv")

# Compute dynamic range
dynamic_range = pd.read_csv("dynamic_range.csv")# same shape (16, 18)
# Number of participants
n_participants = M_level.shape[1]
subjects = ['S22','S29','S38','S40', 'S41','S42','S43','S46','S47','S48','S49R','S50','S52','S53','S54','S55','S56','S57']

# Set up figure
fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Spacing setup
x_spacing = 1.5
x = np.arange(n_participants) * x_spacing
width = 0.5

# --- Plot 1: M_level neurons ---
b1 = axes[0].boxplot(
    [M_neurons.iloc[:, i] for i in range(n_participants)],
    positions=x,
    widths=width,
    patch_artist=True,
    boxprops=dict(facecolor='skyblue', alpha=0.7),
    medianprops=dict(color='black')
)
axes[0].set_ylabel("M_level neurons", color='skyblue')
axes[0].set_title("M_level Neurons Across Participants")

# --- Plot 2: Dynamic Range ---
b2 = axes[1].boxplot(
    [dynamic_range.iloc[:, i] for i in range(n_participants)],
    positions=x,
    widths=width,
    patch_artist=True,
    boxprops=dict(facecolor='lightcoral', alpha=0.7),
    medianprops=dict(color='black')
)
axes[1].set_ylabel("Dynamic Range (M_level - MP) [dB]", color='coral')
axes[1].set_title("Dynamic Range Across Participants")

# --- Common formatting ---
axes[1].set_xticks(x)
axes[1].set_xticklabels(subjects, rotation=45)
axes[1].set_xlabel("Participant")

for ax in axes:
    ax.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 6))

colors = plt.cm.tab20(np.linspace(0, 1, n_participants))  # unique colors per subject

# --- Plot each subject’s data ---
for i, subj in enumerate(subjects):
    x_vals = M_neurons.iloc[:, i].values
    y_vals = dynamic_range.iloc[:, i].values

    # Mask NaNs
    mask = ~np.isnan(x_vals) & ~np.isnan(y_vals)
    x_vals = x_vals[mask]
    y_vals = y_vals[mask]

    # Scatter for each subject
    plt.scatter(
        x_vals, y_vals,
        color=colors[i],
        label=subj,
        s=50, alpha=0.7,
        edgecolor='black', linewidth=0.3
    )

    # Individual regression line
    if len(x_vals) > 1:  # need at least 2 points
        m_i, b_i = np.polyfit(x_vals, y_vals, 1)
        x_fit = np.linspace(min(x_vals), max(x_vals), 50)
        y_fit = m_i * x_fit + b_i
        plt.plot(x_fit, y_fit, color=colors[i], alpha=0.7, linewidth=1.5)

x_all = M_neurons.values.flatten()
y_all = dynamic_range.values.flatten()
mask = ~np.isnan(x_all) & ~np.isnan(y_all)
x_all = x_all[mask]
y_all = y_all[mask]
# --- Group-level regression line (black dashed) ---
r_group, p_group = pearsonr(x_all, y_all)
m, b = np.polyfit(x_all, y_all, 1)
x_fit = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
y_fit = m * x_fit + b
plt.plot(x_fit, y_fit, color='black', linestyle='--', linewidth=3,
         label=f'Group fit (r={r_group:.2f}, p={p_group:.3f})')

# --- Formatting ---
plt.xlabel("M level neurons")
plt.ylabel("Dynamic Range (M_level - MP) dB")
plt.title("M level neurons vs Dynamic Range Across Electrodes and Participants")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# --- Print correlations ---

print(f"\n=== Group-level correlation ===\n r = {r_group:.2f}, p = {p_group:.3f}")

m, b = np.polyfit(x_all, y_all, 1)
x_fit = np.linspace(np.nanmin(x_all), np.nanmax(x_all), 200)
y_fit = m * x_fit + b

plt.plot(x_fit, y_fit, color='black', linestyle='--', linewidth=3, label=f'Group fit: y={m:.2f}x+{b:.2f}')

# --- Formatting ---
plt.xlabel("M level neurons")
plt.ylabel("Dynamic Range (M_level - MP) dB")
plt.title("M level neurons vs Dynamic Range Across Electrodes and Participants")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()