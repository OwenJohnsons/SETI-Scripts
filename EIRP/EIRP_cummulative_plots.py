'''
Code Purpose: 
Author: Owen A. Johnson
Last Major Update: 12/01/2023
'''
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from EIRP_functions import *
import scienceplots
plt.style.use(['science','ieee'])

# --- .csv import --- 
gaia_df = pd.read_csv('gaia_inbeam_table.csv')
print('Amount of Gaia targets loaded:', len(gaia_df))

TESS_df = pd.read_csv('TESS_ext_target_data_observed.csv')
print('Amount of TESS targets loaded:', len(TESS_df))

# --- .txt pointings import --- 

pointings_ra, pointings_dec = np.loadtxt('pointings-10122022.txt', unpack = True)
# gaia_ra, gaia_dec = np.loadtxt('gaia_targets_in_beam_coords_p1.295.txt', unpack = True)
gaia_ra = gaia_df['ra']; gaia_dec = gaia_df['decl']

# --- Plotting --- 

fig, axes = plt.subplots(figsize=(8,10), dpi = 200)
axes.set_aspect(1)

for i in range(0, len(pointings_ra)):
    circle120 = plt.Circle((pointings_ra[i], pointings_dec[i]), 1.295, fill = False, lw = 0.15, color = 'green')
    axes.add_artist(circle120)

plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
plt.scatter(pointings_ra, pointings_dec, s = 0.5, color = 'red', zorder = 0)
# plt.scatter(gaia_ra, gaia_dec, s = 0.1, color = 'black', zorder = 3, alpha = 0.1)
plt.scatter(gaia_df['ra'], gaia_df['decl'], s = 0.05, color = 'black', zorder = 3, alpha = 0.1)
plt.savefig('Gaia_targets_in_beam_pointings.pdf')


plt.title('Gaia Targets within FWHM of LOFAR Beam Pointings')
plt.show()

# --- EIRP calculations ---

# --- Constants for LOFAR Survey --- #
Jansky = 1e-26 # 1 Jansky = 1e-26 Watts/m^2/Hz
obs_time = 15*60 # 15 minutes in seconds
SEFD = 5e3*Jansky # System Equivalent Flux Density for LOFAR-HBA

# --- Transmitter EIRP --- 
EIRP_IP_radar = np.log10(pow(10,17)) # EIRP of a K1-type transmitter
EIRP_LR_radar = np.log10(pow(10,13)) # EIRP of a planetary radar
EIRP_TV_broadcase = np.log10(pow(10,10)) # EIRP of a Aircraft radar

narrowband_gaia = np.log10(obsEIRP(5, SEFD, gaia_df['dist_c'], obs_time, 3, 0.1))
narrowband_TESS = np.log10(obsEIRP(5, SEFD, TESS_df['d'], obs_time, 3, 0.1)) 
print(np.max(narrowband_gaia), np.min(narrowband_gaia))

# --- EIRP Histogram Plot --- 
plt.figure(figsize=(8,6), dpi = 200)
plt.hist(narrowband_gaia, bins = 50, color = 'black', alpha = 0.5, label='Gaia')
plt.hist(narrowband_TESS, bins = 50, color = 'blue', alpha = 0.5, label='TESS')

plt.axvline(EIRP_IP_radar, label="K1-type", linestyle = '--', linewidth = 1)
plt.axvline(EIRP_LR_radar, label="Planetary Radar", linestyle = '-.', linewidth = 1)
plt.axvline(EIRP_TV_broadcase, label="Aircraft Radar", linestyle = ':', linewidth = 1)

plt.legend(loc = 'upper left', frameon = True)
plt.xlabel("EIRP (W/Hz)")
plt.ylabel("Counts")

plt.savefig('EIRP_histogram_plot.pdf')
plt.show()

# --- Cummalative EIRP plot ---

fig, ax = plt.subplots(figsize=(8,6), dpi = 200)

n, bins, patches = plt.hist(narrowband_gaia, bins = 1000, color = 'black', alpha = 0.5, histtype = 'stepfilled', density = True, cumulative = True, label = 'Gaia')
n, bins, patches = plt.hist(narrowband_TESS, bins = 1000, color = 'blue', alpha = 0.5, histtype = 'stepfilled', density = True, cumulative = True, label = 'TESS')

ax.axvline(EIRP_IP_radar, label="K1-type", linestyle = '--', linewidth = 1)
ax.axvline(EIRP_LR_radar, label="Planetary Radar", linestyle = '-.', linewidth = 1)
ax.axvline(EIRP_TV_broadcase, label="Aircraft Radar", linestyle = ':', linewidth = 1)

ax.legend(loc = 'upper left', frameon = True)
ax.set_title('Cumulative EIRP Distribution for Gaia Targets')
ax.set_xlabel('EIRP (W/Hz)'); ax.set_ylabel('Cumulative Probability')

plt.savefig('EIRP_cummulative_plot.pdf')
plt.show()

# --- Distance Normalized Histogram Plot ---

plt.figure(figsize=(8,6), dpi = 200)
plt.hist(gaia_df['dist_c'], bins = 50, color = 'black', alpha = 0.5, label='Gaia', density = True)
plt.hist(TESS_df['d'], bins = 10, color = 'blue', alpha = 0.5, label='TESS', density = True)

plt.axvline(gaia_df['dist_c'].mean(), label = 'Gaia Mean Distance: %s pc' % '{:.1f}'.format(gaia_df['dist_c'].mean()), linestyle = '--', linewidth = 1)
plt.axvline(TESS_df['d'].mean(), label = 'Tess Mean Distance: %s pc' % '{:.1f}'.format(TESS_df['d'].mean()), linestyle = '-.', linewidth = 1)

plt.xlabel('Distance (pc)')
# plt.ylabel('')
plt.legend()