'''
Code Purpose: Query the gaia database for targets in the FWHM of the LOFAR beam pointing. 
Author: Owen A. Johnson
Last Major Update: 19/12/2022
'''

import sqlite3
import time 
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm
import pandas as pd 

connection = sqlite3.connect('gaia_master_database.db')
connection.row_factory = sqlite3.Row
cursor = connection.execute('select * from data') # --- fetching 'data' table from database ---

# --- fetching column names --- 
row = cursor.fetchone()
names = row.keys()

# --- fetching RA and DEC from the database ---

def query_column(column_name):
    start = time.time() # timing the query 
    cursor.execute('SELECT %s FROM data' % column_name) # actual query 
    entries = cursor.fetchall()
    end = time.time()
    print('Time taken to fetch %s from database: %s secs.' % (column_name, "{:.1f}".format(end - start)))
    return entries 

ra_db = query_column('ra')
dec_db = query_column('decl')
sourceid_db = query_column('source_id')

# --- Loading LOFAR beam pointings ---
pointings_ra, pointings_dec = np.loadtxt('/datax/scratch/owenj/pointings-10122022.txt', unpack = True) # Unpacking LOFAR pointings in RA and DEC 
obs_trgts = SkyCoord(pointings_ra, pointings_dec, frame='fk5', unit=(u.hourangle, u.deg)) # LOFAR beam pointing RA and DEC

# --- Looping pointings through Gaia database --- 

FWHM_tois = [] # Empty array to store TOIs in FWHM of LOFAR beam pointing

for i in tqdm(range(0, len(obs_trgts))): 
    pointing = obs_trgts[i]
    for j in tqdm(range(0, len(ra_db)), leave = False):
        gaia_trgt = SkyCoord(ra_db[j][0], dec_db[j][0], frame='fk5', unit=(u.hourangle, u.deg))
        sep = pointing.separation(gaia_trgt)

        if sep < 2.59*u.deg: # FWHM of LOFAR beam from Van Haarlem et al. 2013
            print('Target found in beam pointing: %s' % pointing)
            print('Target found in database: %s' % sourceid_db[j][0])
            print('Separation: %s' % sep.degree)
            FWHM_tois.append((sourceid_db[j][0], sep.degree, gaia_trgt.ra.degree, gaia_trgt.dec.degree, pointing.ra.degree, pointing.dec.degree))

# --- Saving TOIs to a dataframe ---
FWHM_tois_df = pd.DataFrame(FWHM_tois, columns = ['Gaia Source ID', 'Separation (Deg)', 'Source RA (deg)', 'Source DEC (deg)', 'Pointing RA (deg)', 'Pointing DEC (deg)'])
FWHM_tois_df.to_csv('FWHM_tois.csv', index = False)