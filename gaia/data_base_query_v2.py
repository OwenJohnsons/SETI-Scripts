'''
Code Purpose: Query the gaia database for targets in the FWHM of the LOFAR beam pointing. New version using vectorisation to speed up the code as suggested by Evan. 
Author: Owen A. Johnson
Last Major Update: 09/01/2023
'''

import sqlite3
import time 
import numpy as np
np.set_printoptions(suppress=True)
from tqdm import tqdm
import pandas as pd 

script_runtime_start = time.time()

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

ra_db = query_column('ra')
dec_db = query_column('decl')
# sourceid_db = query_column('source_id')
targets_vec = np.vstack((ra_db, dec_db)).T

# --- Loading LOFAR beam pointings ---
pointings_ra, pointings_dec = np.loadtxt('/datax/scratch/owenj/pointings-10122022.txt', unpack = True) # Unpacking LOFAR pointings in RA and DEC 
pointings_vec = np.vstack((pointings_ra, pointings_dec)).T

# --- Finding targets in beam --- 

def FWHM_seperation(points, beam_pointing, seperation_limit):
    points = np.array(points)
    beam_pointing = np.array(beam_pointing)
    seperation = np.linalg.norm(points - beam_pointing, axis=1)
    return points[seperation <= seperation_limit]

sep_limit = 2.59 # Taken as the beam FWHM from Van Haarlem et al. 2013

target_count = 0
for pointing in tqdm(pointings_vec):
    in_beam_targets = FWHM_seperation(targets_vec, pointing, sep_limit)
    target_count += len(in_beam_targets)

print('Number of Gaia targets found in beam: %s' % target_count)

script_runtime_end = time.time()
print('Time taken for script to run %s secs.' % ("{:.1f}".format(script_runtime_start - script_runtime_end)))