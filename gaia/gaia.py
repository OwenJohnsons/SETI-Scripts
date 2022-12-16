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

# - Looping over LOFAR beam pointings to find Gaia targets within the FWHM of the beam -

in_beam_targets = np.array([])

for i in tqdm(range(0, len(obs_trgts))):
    pointing = obs_trgts[i]
    idx = 0 
    for j in range(0, len(gaia_df)):
        idx += 1 
        gaia_trgt = SkyCoord(gaia_df['ra'][j], gaia_df['dec'][j], frame='fk5', unit=(u.hourangle, u.deg))
        sep = pointing.separation(gaia_trgt)

        if sep < 2.59*u.deg: # FWHM of LOFAR beam from Van Haarlem et al. 2013
            in_beam_targets = np.append(in_beam_targets, (gaia_df['ra'][j], gaia_df['dec'][j]))

# - Saving targets in beam to a text file -
np.savetxt('in_beam_targets.txt', in_beam_targets, fmt='%s')

gaia_df.visualize(filename='gaia_df_task_graph.pdf') 