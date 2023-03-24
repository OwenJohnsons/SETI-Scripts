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
import pickle
from astropy import units as u
from astropy.coordinates import SkyCoord

def calc_DishArea(d):
    """ Compute dish area
    d = dish diameter
    """
    return np.pi * (d/2)**2
    
# --- .csv import --- 
gaia_df = pd.read_csv('csv/gaia_merged.csv').sort_values(by=['beam_ra'])
print('Amount of Gaia targets loaded:', len(gaia_df))

TESS_df = pd.read_csv('csv/TESS_ext_target_data_observed.csv').sort_values(by=['ra'])
TESS_ids = TESS_df['TIC_ID'].values
print('Amount of TESS targets loaded:', len(TESS_df))

total_targets = len(gaia_df) + len(TESS_df)
print('Total amount of targets:', total_targets)

# --- Matching Beam Number to TIC ID --- 
ided_array = []
for i in range(0, len(TESS_ids)):
    replace_ids = np.where(gaia_df['beam'] == i)
    fill = np.full(len(replace_ids[0]), TESS_ids[i])
    ided_array.append(fill)
gaia_df['beam'] = np.concatenate(ided_array)
# gaia_df.to_csv('csv/gaia_merged_t.csv', index=False)

# --- .txt pointings import --- 
pointings_ra, pointings_dec = np.loadtxt('txt/pointings-10122022.txt', unpack = True)
# gaia_ra, gaia_dec = np.loadtxt('gaia_targets_in_beam_coords_p1.295.txt', unpack = True)
gaia_ra = gaia_df['ra_x']; gaia_dec = gaia_df['dec']

# --- Plotting --- 
fig, axes = plt.subplots(figsize=(8,10), dpi = 200)
axes.set_aspect(1)

for i in range(0, len(pointings_ra)):
    circle120 = plt.Circle((pointings_ra[i], pointings_dec[i]), 1.295, fill = False, lw = 0.15, color = 'green')
    axes.add_artist(circle120)

plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
plt.scatter(pointings_ra, pointings_dec, s = 0.5, color = 'red', zorder = 0)
# plt.scatter(gaia_ra, gaia_dec, s = 0.1, color = 'black', zorder = 3, alpha = 0.1)
plt.scatter(gaia_df['ra_x'], gaia_df['dec'], s = 0.05, color = 'black', zorder = 3, alpha = 0.1)
plt.savefig('plots/Gaia_targets_in_beam_pointings.pdf')
plt.title('Gaia Targets within FWHM of LOFAR Beam Pointings')
plt.show()

# --- EIRP calculations ---

# --- Constants for LOFAR Survey --- #
Jansky = 1e-26 # 1 Jansky = 1e-26 Watts/m^2/Hz
obs_time = 15*60 # 15 minutes in seconds
SEFD = 5e3*Jansky # System Equivalent Flux Density for LOFAR-HBA

# # --- Transmitter EIRP --- 
EIRP_IP_radar = np.log10(pow(10,17)) # EIRP of a K1-type transmitter
EIRP_LR_radar = np.log10(pow(10,13)) # EIRP of a planetary radar
EIRP_TV_broadcase = np.log10(pow(10,10)) # EIRP of a Aircraft radar
EIRP_EARTH = np.log10(0.7*pow(10, 16)) # EIRP of the Earth

narrowband_gaia = np.log10(obsEIRP(5, SEFD, gaia_df['dist_c'], obs_time, 3, 0.1))
narrowband_TESS = np.log10(obsEIRP(5, SEFD, TESS_df['d'], obs_time, 3, 0.1)) 
print(np.max(narrowband_gaia), np.min(narrowband_gaia))

# # --- EIRP Histogram Plot --- 
plt.figure(figsize=(8,6), dpi = 200)
plt.hist(narrowband_gaia, bins = 50, color = 'black', alpha = 0.5, label='Gaia')
plt.hist(narrowband_TESS, bins = 50, color = 'blue', alpha = 0.5, label='TESS')

plt.axvline(EIRP_IP_radar, label="K1-type", linestyle = '--', linewidth = 1)
plt.axvline(EIRP_LR_radar, label="Planetary Radar", linestyle = '-.', linewidth = 1)
plt.axvline(EIRP_TV_broadcase, label="Aircraft Radar", linestyle = ':', linewidth = 1)

plt.legend(loc = 'upper left', frameon = True)
plt.xlabel("EIRP (W/Hz)")
plt.ylabel("Counts")

plt.savefig('plots/EIRP_histogram_plot.pdf')
plt.show()

# # --- Cummalative EIRP plot ---

fig, ax = plt.subplots(figsize=(8,3.5), dpi = 200)

n, bins, patches = plt.hist(narrowband_gaia, bins = 1000, color = 'black', alpha = 0.5, histtype = 'stepfilled', density = True, cumulative = True, label = 'Gaia')
n, bins, patches = plt.hist(narrowband_TESS, bins = 1000, color = 'blue', alpha = 0.5, histtype = 'stepfilled', density = True, cumulative = True, label = 'TESS')

ax.axvline(EIRP_IP_radar, label="K1-type", linestyle = '--', linewidth = 1)
ax.axvline(EIRP_EARTH, label="Earth-type", linewidth = 1)
ax.axvline(EIRP_LR_radar, label="Planetary Radar", linestyle = '-.', linewidth = 1)
ax.axvline(EIRP_TV_broadcase, label="Aircraft Radar", linestyle = ':', linewidth = 1)

ax.legend(loc = 'upper left', frameon = True)
# ax.set_title('Cumulative EIRP Distribution for Gaia Targets')
ax.set_xlabel('EIRP (W/Hz)'); ax.set_ylabel('Cumulative Probability')

plt.savefig('plots/EIRP_cummulative_plot.pdf')
plt.show()

# # --- Distance Normalized Histogram Plot ---
plt.figure(figsize=(8,6), dpi = 200)
plt.hist(gaia_df['dist_c'], bins = 50, color = 'black', alpha = 0.5, label='Gaia', density = True)
plt.hist(TESS_df['d'], bins = 10, color = 'blue', alpha = 0.5, label='TESS', density = True)

plt.axvline(gaia_df['dist_c'].mean(), label = 'Gaia Mean Distance: %s pc' % '{:.1f}'.format(gaia_df['dist_c'].mean()), linestyle = '--', linewidth = 1)
plt.axvline(TESS_df['d'].mean(), label = 'Tess Mean Distance: %s pc' % '{:.1f}'.format(TESS_df['d'].mean()), linestyle = '-.', linewidth = 1)

plt.xlabel('Distance (pc)')
# plt.ylabel('')
plt.legend()
plt.show()

#%%
# --- Tsys Histograms --- 
Tsys_df = pd.read_csv('Tsys/Tsys.csv')
# print(Tsys_df.head())
freqs = np.arange(110, 200, 10)

    
# --- Function Plots --- 

plt.figure(figsize=(10,6), dpi = 200)
for id in Tsys_df['TICid'].unique():
    Tsys_indv = Tsys_df[Tsys_df['TICid'] == id]

    freqs = Tsys_indv['freq']
    tt = Tsys_indv['Tsys']
    
    plt.scatter(freqs, tt, label = 'TIC %s' % id, s = 2)
    plt.plot(freqs, tt, linestyle = '--', linewidth = 1, alpha = 0.1)
    
plt.xlabel('Frequency (MHz)')
plt.ylabel('Tsys (K)')
# plt.legend()
plt.savefig('plots/TsysvsFreq_plot.png')
plt.show()

lofar_dish_area = 1677.6 #calc_DishArea(100)

plt.figure(figsize=(10,6), dpi = 200)
for id in Tsys_df['TICid'].unique():
    Tsys_indv = Tsys_df[Tsys_df['TICid'] == id]
    SEFD = calc_SEFD(lofar_dish_area, Tsys_indv['Tsys'], eff=1.0)
    freqs = Tsys_indv['freq']
    
    plt.scatter(freqs, SEFD, label = 'TIC %s' % id, s = 2)
    plt.plot(freqs, SEFD, linestyle = '--', linewidth = 1, alpha = 0.1)
    
plt.xlabel('Frequency (MHz)')
plt.ylabel('SEFD (Jy)')
# plt.legend()
plt.savefig('plots/SEFDvsFreq_plot.png')
plt.show()


# --- Specific Tsys Histogram Plots (3 x 3) ---

overall_Tsys = []
for freq in freqs:
    Tsys_indv = Tsys_df[Tsys_df['freq'] == freq]
    frequency_values = np.array([])
    for id in Tsys_df['TICid'].unique():
        target_numbers = len(np.where(gaia_df['beam'] == id)[0])
        Tsys_value = Tsys_indv[Tsys_indv['TICid'] == id]['Tsys'].values[0]
        list_values = [Tsys_value] * target_numbers
        # print(len(list_values), Tsys_value, np.mean(np.array(list_values)))
        frequency_values = np.append(frequency_values, list_values)
    overall_Tsys.append(frequency_values)
    

# --- Specific Tsys Plot: Verticle Plots --- 
#%%
fig, ax = plt.subplots(3, 3, sharex='col', sharey='row', figsize=(10,3), dpi = 500)
freqs = np.arange(110, 200, 10)

for i in range(3):
    for j in range(3):
        ax[i, j].hist(overall_Tsys[i*3+j], bins = 20, color = 'black', alpha = 0.5, density = False, facecolor = 'none', edgecolor = 'black', linewidth = 1, label='%s MHz' % '{:0f}'.format(freqs[i*3+j]))
        # ax[i, j].title.set_text('Frequency: %s MHz' % '{:.0f}'.format(freqs[i*3+j]))
        # ax[i,j].hist(test['Tsys'], bins = 50, color = 'black', alpha = 0.5, label='Gaia', density = True)
        # - plot text in the top right corner of the plot - 
        ax[i, j].text(0.96, 0.85, '%s MHz' % '{:.0f}'.format(freqs[i*3+j]), horizontalalignment='right', verticalalignment='top', transform=ax[i, j].transAxes)

for j in range(3):
    ax[2, j].set_xlabel('$T_{sys}$ (K)')
    ax[j, 0].set_ylabel('Counts')
    
plt.savefig('plots/Tsys_hist_plot.png', bbox_inches='tight')

fig, ax = plt.subplots(3, 3, sharex='col', sharey='row', figsize=(10,3), dpi = 500)

for i in range(3):
    for j in range(3):
        ax[i, j].hist(calc_SEFD(lofar_dish_area, overall_Tsys[i*3+j], eff = 1.0), bins = 20, color = 'black', alpha = 0.5, density = False, facecolor = 'none', edgecolor = 'black', linewidth = 1, label='%s MHz' % '{:0f}'.format(freqs[i*3+j]))
        # ax[i, j].title.set_text('Frequency: %s MHz' % '{:.0f}'.format(freqs[i*3+j]))
        # ax[i,j].hist(test['Tsys'], bins = 50, color = 'black', alpha = 0.5, label='Gaia', density = True)
        # - plot text in the top right corner of the plot - 
        ax[i, j].text(0.96, 0.85, '%s MHz' % '{:.0f}'.format(freqs[i*3+j]), horizontalalignment='right', verticalalignment='top', transform=ax[i, j].transAxes)

for j in range(3):
    ax[2, j].set_xlabel('SEFD (Jy)')
    ax[j, 0].set_ylabel('Counts')
    
plt.savefig('plots/SEFD_hist_plot.png', bbox_inches='tight')


# --- Sensitivity Plot ---
#%%
fig, ax = plt.subplots(3, 3, sharex='col', sharey='row', figsize=(10,3), dpi = 500)
for i in range(3):
    for j in range(3):
        SEFD = calc_SEFD(lofar_dish_area, overall_Tsys[i*3+j], eff = 1.0)
        narrowband = np.log10(obsEIRP(5, SEFD*Jansky, gaia_df['dist_c'], obs_time, 3, 0.1))
        ax[i, j].hist(narrowband, bins = 100, color = 'black', alpha = 0.5, facecolor = 'none', edgecolor = 'black', linewidth = 1, histtype = 'stepfilled', density = True, cumulative = True)
        ax[i, j].text(0.25, 0.85, '%s MHz' % '{:.0f}'.format(freqs[i*3+j]), horizontalalignment='right', verticalalignment='top', transform=ax[i, j].transAxes)
        ax[i, j].axvline(EIRP_IP_radar, label="K1-type", linestyle = '--', linewidth = 1)
        ax[i, j].axvline(EIRP_LR_radar, label="Planetary Radar", linestyle = '-.', linewidth = 1)
        ax[i, j].axvline(EIRP_TV_broadcase, label="Aircraft Radar", linestyle = ':', linewidth = 1)

for j in range(3):
    ax[2, j].set_xlabel('EIRP (W/Hz)')
    # ax[j, 0].set_ylabel('Cumulative Probability')
    
plt.savefig('plots/EIRP_hist_plot.png', bbox_inches='tight')


# --- Galactic Coordinates vs. Tsys ---
#%%
galactic_coords_tess = SkyCoord(ra = TESS_df['ra']*u.degree, dec = TESS_df['dec']*u.degree, frame = 'icrs').galactic

Tsys_indv = Tsys_df[Tsys_df['TICid'] == id]

plt.figure(figsize=(8, 4), dpi = 500)
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.scatter(galactic_coords_tess.l.wrap_at(180*u.degree).radian, galactic_coords_tess.b.radian, c = Tsys_df[Tsys_df['freq'] == 190]['Tsys'].values, cmap = 'viridis', s = 5)
plt.title('Frequency: %s MHz' % 190)
plt.axhspan(np.deg2rad(-5), np.deg2rad(5), alpha=0.2, color='grey', label = 'Galactic Center')

plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.legend()
values = Tsys_df[Tsys_df['freq'] == 190]['Tsys'].values
# xticks = np.arange(np.round(np.min(values), 0), np.round(np.max(values), 15, 0))
# plt.colorbar(orientation="horizontal").ax.set_xticklabels(xticks, rotation=45)
plt.colorbar(label = 'Tsys (K)', orientation="horizontal")
plt.savefig('plots/Tsys_galactic_coords_190.png', bbox_inches='tight', dpi = 500)
plt.show()
