'''
Code Purpose: 
Author: Owen A. Johnson
Last Major Update: 01/02/2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.stats as stats
import scienceplots
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit
plt.style.use(['science','ieee'])


def normpdf(x, mean, sd):
    var = float(sd)**2
    denom = (2*np.pi*var)**.5
    num = np.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom


fwhm = 2.59 #deg
mu = 0
sigma = fwhm/2.355
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)

scale_factor = 1/normpdf(0, mu, sigma)

# --- Fitting Exp function --- 

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

smooth_x = np.linspace(0, 1.295, 1000)
popt, pcov = curve_fit(func, xdata, ydata)

plt.plot(x, stats.norm.pdf(x, mu, sigma)*scale_factor)
plt.axvspan(mu - fwhm/2, mu + fwhm/2, alpha=0.2, label = 'FWHM')
plt.axvline(0, color='red', label = 'Pointing')

plt.title('$ \mathrm {FWHM} =2{\sqrt {2\ln 2}}\;\sigma \\approx 2.355\;\sigma$')
plt.legend()
plt.show()

# --- Plot for Paper --- 

OB_angle = np.arange(0.1, 1, 0.1)*1.295
angles = np.linspace(0, 1.295, 1000); sensitivity = stats.norm.pdf(angles, mu, sigma)*scale_factor

files = sorted(glob.glob('data/gaia_targets_in_beam_coords_p*.txt'))

target_no = []
for file in files:
    target_no.append(np.loadtxt(file, unpack = True)[0].size)

fig, ax1 = plt.subplots(figsize = (4, 4))
plt.ticklabel_format(axis='y', style='sci', scilimits=(4,4))

ax2 = ax1.twinx()
ax2.scatter(OB_angle[0:len(files)], target_no, s = 10, color = 'black', facecolor = 'none', edgecolors='black')
ax2.scatter(1.295, 96498, marker = 'x', color = 'red', s = 10, label = 'FWHM')
ax2.set_ylabel('Number of Gaia Targets', rotation=-90, labelpad= 10)
ax2.set_yscale('log')

ax1.set_xlabel('Off-Beam Angle, $\\theta$ [deg$^\circ$ ]')
ax1.set_ylabel('Beam Sensitivity')
ax1.legend(loc = 'center left'); ax2.legend(loc = 'center left')
ax1.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=None))
ax1.plot(angles, sensitivity*100, 'b-', label = 'Sensitivity Curve', alpha = 0.5)

plt.savefig('off-beam-sensitivity-plots.pdf', bbox_inches = 'tight')
plt.show()
# %%

# %%
