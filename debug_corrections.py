# script to check on algorithm for correcting electrode distance fits
import numpy as np
import matplotlib.pyplot as plt


def mean_error(arg, sl, intercept):
    return (arg * sl) + intercept


# takes as argument the fit distance and the average slope and intercept for errors
# from a population of subjects and returns a corrected fit distance, constrained to between 0.05 and 1.95 mm
def correct_error(f_d, sl, inter):
    err = (f_d * sl) + inter
    corrected_dist = f_d - err
    return corrected_dist, err


# make an array of 10 equally spaced numbers, corresponding to ct_dist
ct_dist = np.arange(10)  # zero through 9, just to get a spread

# assume some linear relationship
slope = 0.15
intercept = 0.3

fit_dist = ct_dist * slope + intercept  # set fake fit distances

# add noise
noise = np.random.rand(10)
fit_dist += noise

err_dist = fit_dist - ct_dist  # calculate errors

fig, ax = plt.subplots(3, 2, figsize=(8, 10))
fig.tight_layout(pad=3)
ax[0, 0].plot(ct_dist, fit_dist, '.')  # Plot fit dist v. ct dist
ax[0, 0].set_xlabel('ct dist')
ax[0, 0].set_ylabel('fit_dist')

ax[0, 1].plot(fit_dist, err_dist, '.')  # Plot error v. fit dist, as in paper
ax[0, 1].set_xlabel('ct dist')
ax[0, 1].set_ylabel('err_dist')

[err_sl, err_int] = np.polyfit(fit_dist, err_dist, 1)  # linear fit slope and intercept
ax[0, 1].plot((fit_dist[0], fit_dist[-1]), (fit_dist[0] * err_sl + err_int, fit_dist[-1] * err_sl + err_int))
ax[0, 1].set_xlabel('fit dist')

dist_corr, err = correct_error(fit_dist, err_sl, err_int)  # Call function to calculate correction

# err = fit_dist*err_sl + err_int  # leftover code from before I created the function

ax[1, 0].plot(fit_dist, err, '.')
ax[1, 0].set_xlabel('fit dist')
ax[1, 0].set_ylabel('err')

# dist_corr = fit_dist - err  # leftover code from before I created the function
ax[1, 1].plot(ct_dist, dist_corr, '.')
ax[1, 1].set_xlabel('ct dist')
ax[1, 1].set_ylabel('dist_corr')

err_after_corr = dist_corr - ct_dist
# error after correction
ax[2, 0].plot(ct_dist, err_after_corr, '.')
ax[2, 0].set_xlabel('ct dist')
ax[2, 0].set_ylabel('err_after_corr')

plt.show()
