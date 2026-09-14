import matplotlib.pyplot as plt
import numpy as np

elecs = np.arange(16) + 1
thresh_db = np.array([33.16878491, 33.16878491, 33.16878491, 32.74499893, 32.74499893, 32.29947104, 32.74499893,
            32.29947104, 33.57285263, 33.57285263, 33.57285263, 33.16878491, 32.74499893, 32.74499893,
            32.74499893, 32.74499893])
mcl_db = np.array([53.04380593, 53.04380593, 53.04380593, 52.78832016, 52.95947624, 52.39042359, 52.52509124, 52.4355452,
            53.76804938, 53.41347404, 53.57285263, 53.08566572, 52.52509124, 52.74499893, 52.74499893, 52.74499893])
rpos = np.array([0.87769811, 0.79012327, 0.39010671, -0.02793052, -0.22808854, -0.25158784, -0.1396573, 0.04036908,
            0.21766321, 0.52362401, 0.85145473, 0.82622878, 0.74981766, 0.50357691, 0.32183185, -0.03709997])
dist = 1-rpos
surv_vals = [0.50, 0.32, 0.18, 0.16, 0.08, 0.00, 0.05, 0.17, 0.24, 0.11, 0.16, 0.17, 0.34,
            0.16, 0.24, 0.04]

# fig, axs = plt.subplots(nrows=2, ncols=1)
fig, axs = plt.subplots(3, 1, figsize=(6, 10))
# fig.suptitle('Subject 2596')

axs[0].plot(elecs, dist, 'ob', label='distance (mm)')
axs[0].set_ylabel('distance (mm)')

axs[1].plot(elecs, surv_vals, 'or', label='survival')
axs[1].set_ylabel('survival')



axs[2].plot(elecs, thresh_db, 'or', label='Threshold')
axs[2].plot(elecs, mcl_db, 'ob', label='MCL')
axs[2].set_xlabel('Electrode')
axs[2].set_ylabel('Stim intensity (dB)')
axs[2].set_ylim([30, 60])
axs[2].legend()


plt.show()
