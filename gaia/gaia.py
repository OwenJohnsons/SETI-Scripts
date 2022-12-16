"""
Code Purpose: Using DASK to see targets in the Gaia database that fall within the FWHM of the LOFAR beam pointing at a given RA and DEC. This code contains two hardcoded file paths that need to be changed to the correct paths on the user's machine. In this case the paths are specified on the Breakthrough Listen (blpc1) node. The database is approximately 13 GB in size thus use caution when conventionally loading the data base using something like pd.read_csv().
Author: Owen A. Johnson 
Last Major Update: 14/12/2022
"""

import time 
import numpy as np
import pandas as pd 
from dask import dataframe as dd 
import astropy 
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u

# - Hardcoded file paths - 
gaia_database_path = '/datax/scratch/vgajjar/master_gaia_database.csv'
pointings_ra, pointings_dec = np.loadtxt('/datax/scratch/owenj/pointings-10122022.txt', unpack = True) # Unpacking LOFAR pointings in RA and DEC 

# - Loading beam pointing data -
obs_trgts = SkyCoord(pointings_ra, pointings_dec, frame='fk5', unit=(u.hourangle, u.deg)) # LOFAR beam pointing RA and DEC

# - Loading Gaia database using DASK -
start = time.time()
gaia_df = dd.read_csv(gaia_database_path) # 1 million rows at a time
end = time.time()
print('Time to load Gaia database: %s secs' % (end - start)) 
print(gaia_df.head())

# - Empty Dataframe to store targets in beam with Gaia headers -
in_beam_targets_df = pd.DataFrame(columns = ['ra', 'decl', 'parallax', 'pmra', 'pmdec', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'radial_velocity', 'radial_velocity_error', 'l', 'b', 'ecl_lon', 'ecl_lat', 'teff_val', 'a_g_val', 'e_bp_min_rp_val', 'radius_val', 'lum_val', 'datalink_url'], dtype = 'object')

# - Looping over LOFAR beam pointings to find Gaia targets within the FWHM of the beam -

in_beam_targets = np.array([])

for i in tqdm(range(0, len(obs_trgts))):
    pointing = obs_trgts[i]
    for idx,row in gaia_df.iterrows(): # Iterating each pointing over the gaia database
        gaia_trgt = SkyCoord(row['ra'], row['decl'], frame='fk5', unit=(u.hourangle, u.deg))
        sep = pointing.separation(gaia_trgt)

        if sep < 2.59*u.deg: # FWHM of LOFAR beam from Van Haarlem et al. 2013
            in_beam_targets_df.append(row)
 
# - Saving targets in beam to a text file -
in_beam_targets_df.to_csv('in_beam_targets.csv')

gaia_df.visualize(filename='gaia_df_task_graph.pdf') 