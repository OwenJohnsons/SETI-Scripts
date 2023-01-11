'''
Code Purpose: Query the gaia database for targets in the FWHM of the LOFAR beam pointing. New version using vectorisation to speed up the code as suggested by Evan. 
Author: Owen A. Johnson
Last Major Update: 09/01/2023
'''

import argparse
import sqlite3
import time 
import numpy as np
from tqdm import tqdm
import pandas as pd 
import matplotlib.pyplot as plt

script_runtime_start = time.time()

# --- Arguments --- 

parser = argparse.ArgumentParser(description='Query the gaia database for targets in the FWHM of the LOFAR beam pointing.')
parser.add_argument('-p', '--percentage', type=float, help='Percentage of LOFAR FWHM beam to be considered suitable for targets, default = 1.'); parser.set_defaults(percentage=1)
percentage = parser.parse_args().percentage

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
    # cursor.execute("SELECT %s FROM data LIMIT 100000" % column_name) # limits the query to 'n' entries; for testing purposes. 
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

sep_limit = 1.295*percentage # Taken as the beam FWHM from Van Haarlem et al. 2013
total_targets = []; target_count = 0

for pointing in tqdm(pointings_vec):
    in_beam_targets = FWHM_seperation(targets_vec, pointing, sep_limit)
    if len(in_beam_targets) != 0: 
        total_targets.append(in_beam_targets) # - appending to a (2, n) array shape 
        target_count += len(in_beam_targets)

# --- Plotting the results ---

fig, axes = plt.subplots(figsize=(8,10), dpi = 200)
axes.set_aspect(1)

for pointing in total_targets:
    plt.scatter(pointing[:,0], pointing[:,1], s=0.1, c='k', alpha=0.001)

for i in range(0, len(pointings_vec)):
    circle120 = plt.Circle((pointings_vec[i]), 2.59, fill = False, lw = 0.15, color = 'green')
    axes.add_artist(circle120)

plt.scatter(pointings_ra, pointings_dec, s=0.2, c='r', label = 'Beam Pointings')
plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
plt.legend()
plt.title('Gaia Targets within FWHM of LOFAR Beam Pointings')
plt.savefig('gaia_targets_in_beam.png', dpi=300, bbox_inches='tight')

# --- Printing results to terminal & exporting ---

print('Number of Gaia targets found in beam: %s' % target_count)

# --- Exporting the results to a text file ---

total_targets = np.array(total_targets, dtype=object)

flatten_array = []
for i in range(0, len(total_targets)):
    for j in range(0, len(total_targets[i])):
        flatten_array.append(total_targets[i][j])

flatten_array = np.array(flatten_array, dtype=object)
np.savetxt(('gaia_targets_in_beam_coords_p%s.txt' % percentage), flatten_array, fmt='%s')

script_runtime_end = time.time()
print('Time taken for script to run %s secs.' % ("{:.1f}".format(script_runtime_end - script_runtime_start)))