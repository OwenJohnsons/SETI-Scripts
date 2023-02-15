import argparse
import sqlite3
import time 
import numpy as np
from tqdm import tqdm
import pandas as pd 
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

script_runtime_start = time.time()

# --- Connecting to database ---
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
    # cursor.execute("SELECT %s FROM data LIMIT 1000000" % column_name) # limits the query to 'n' entries; for testing purposes. 
    entries = cursor.fetchall()
    end = time.time()
    print('Time taken to fetch %s from database: %s secs.' % (column_name, "{:.1f}".format(end - start)))
    return np.array(entries).flatten()

sourceid_db = query_column('source_id')
ra_db = query_column('ra'); dec_db = query_column('decl')
targets_vec = np.vstack((ra_db, dec_db)).T
print('Total number of targets loaded from the Gaia database: %s' % len(targets_vec))

# --- Loading LOFAR beam pointings ---
OBS_df = pd.read_csv('/datax/scratch/owenj/SETI-Scripts/TESS/TESS_ext_target_data_observed.csv') # Loading observed targets
pointings_ra = OBS_df['ra']; pointings_dec = OBS_df['dec']
pointings_vec = np.vstack((pointings_ra, pointings_dec)).T

# --- Using SkyCoord to find Seperation between Pointings and Targets ---
pointing_coords = SkyCoord(ra=pointings_vec[:,0]*u.degree, dec=pointings_vec[:,1]*u.degree)
target_coords = SkyCoord(ra=targets_vec[:,0]*u.degree, dec=targets_vec[:,1]*u.degree)
in_beam_targets = pointing_coords.search_around_sky(target_coords, 1.295*u.degree)
sources = sourceid_db[in_beam_targets[0]]
in_beam_targets = target_coords[in_beam_targets[0]]
print('Number of targets found in beam: %s' % len(np.unique(sources)))
# print(in_beam_targets.ra.degree, in_beam_targets.dec.degree)

# --- Plotting the results ---
plt.scatter(pointings_vec[:,0], pointings_vec[:,1], s=1, c='r', label='Pointings')
plt.scatter(in_beam_targets.ra.degree, in_beam_targets.dec.degree, s=1, c='b', label='Targets', alpha = 0.0001)
plt.savefig('plots/v3_gaia_query.png', dpi=300)

# --- Saving ids to file ---
np.savetxt('txt/gaia_ids.txt', np.unique(sources), fmt='%s')

script_runtime_end = time.time()
print('Time taken for script to run %s secs.' % ("{:.1f}".format(script_runtime_end - script_runtime_start)))