'''
Code Purpose: 
Author: Owen A. Johnson
Last Major Update: 12/01/2023
'''
#%%

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import argparse 
from tqdm import tqdm
from astropy import units as u
from astropy.coordinates import SkyCoord 

# --- Arguments --- 
# parser = argparse.ArgumentParser(description='Query the gaia database for targets in the FWHM of the LOFAR beam pointing.')
# parser.add_argument('-p', '--percentage', type=float, help='Percentage of LOFAR FWHM beam to be considered suitable for targets, default = 1.'); parser.set_defaults(percentage=1)
# percentage = parser.parse_args().percentage

# --- Finding targets in beam --- 
def FWHM_seperation(points, beam_pointing, seperation_limit):
    points = np.array(points)
    beam_pointing = np.array(beam_pointing)
    seperation = np.linalg.norm(points - beam_pointing, axis=1)
    return points[seperation <= seperation_limit]

# --- Loading LOFAR beam pointings ---
OBS_df = pd.read_csv('TESS_ext_target_data_observed.csv')
obs_ra = OBS_df['ra']; obs_dec = OBS_df['dec']
pointings_vec = np.vstack((obs_ra, obs_dec)).T

# --- Loading TESS database --- 
TESS_df = pd.read_csv('TESS_ext_target_data.csv')
# remove the matching IDs from the TESS_df
TESS_df = TESS_df[~TESS_df['TIC_ID'].isin(OBS_df['TIC_ID'])]
print(len(TESS_df))
TESS_ra = TESS_df['ra']; TESS_dec = TESS_df['dec']
targets_vec = np.vstack((TESS_ra, TESS_dec)).T

print(targets_vec)

sep_limit = 1.295*1 # Taken as the beam FWHM from Van Haarlem et al. 2013
total_targets = []; target_count = 0

for pointing in tqdm(pointings_vec):
    in_beam_targets = FWHM_seperation(targets_vec, pointing, sep_limit)
    if len(in_beam_targets) != 0: 
        total_targets.append(in_beam_targets) # - appending to a (2, n) array shape 
        target_count += len(in_beam_targets)
        
        
print('Number of TESS targets in database: ', len(targets_vec))
print('Number of TESS targets in beam: ', target_count)

# --- Plotting the results ---

fig, axes = plt.subplots(figsize=(10,10), dpi = 200)
axes.set_aspect(1)

for pointing in total_targets:
    plt.scatter(pointing[:,0], pointing[:,1], s=0.1, c='k', alpha=0.5)

for i in range(0, len(pointings_vec)):
    circle120 = plt.Circle((pointings_vec[i]), 2.59, fill = False, lw = 0.15, color = 'green')
    axes.add_artist(circle120)

plt.scatter(obs_ra, obs_dec, s=0.2, c='r', label = 'Beam Pointings')

plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
plt.legend()
plt.title('TESS Targets within FWHM of LOFAR Beam Pointings')
plt.savefig('TESS_targets_in_beam.png', dpi=300, bbox_inches='tight')
plt.show()