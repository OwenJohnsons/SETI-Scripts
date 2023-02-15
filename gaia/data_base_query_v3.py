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
print('\nTotal number of targets loaded from the Gaia database: %s' % len(targets_vec))

# --- Loading LOFAR beam pointings ---
OBS_df = pd.read_csv('/datax/scratch/owenj/SETI-Scripts/TESS/TESS_ext_target_data_observed.csv') # Loading observed targets
pointings_ra = OBS_df['ra']; pointings_dec = OBS_df['dec']
pointings_vec = np.vstack((pointings_ra, pointings_dec)).T

# --- Using SkyCoord to find Seperation between Pointings and Targets ---
pointing_coords = SkyCoord(ra=pointings_vec[:,0]*u.degree, dec=pointings_vec[:,1]*u.degree)
target_coords = SkyCoord(ra=targets_vec[:,0]*u.degree, dec=targets_vec[:,1]*u.degree)
in_beam_targets = pointing_coords.search_around_sky(target_coords, 1.295*u.degree)
'''
---
search_around_sky() 
--- 
Key     Function 
0.      indexes of gaia objects around the given pointings
1.      index of pointings that each gaia object is around
2.      On sky distance between the gaia object and the pointing in degrees
3.      3d distance between the gaia object and the pointing in degrees
'''
sources = sourceid_db[in_beam_targets[0]]

print('\n---\nNumber of unique targets found in beam: %s' % len(np.unique(sources)))
print('Number of targets that are seen in multiple beams: %s' % (len(sources) - len(np.unique(sources))))
print('Average distance between targets and pointings: %s deg' % np.round(np.mean(in_beam_targets[2].degree), 3))
print('Max distance between target and a pointing: %s deg' % np.round(np.max(in_beam_targets[2].degree), 3))
print('Pointing with the most Gaia targets around it: %s \n---\n ' % np.argmax(np.bincount(in_beam_targets[1])))


# --- Plotting the results ---
plot_targets = target_coords[in_beam_targets[0]]

fig, axes = plt.subplots(figsize=(8,10), dpi = 200)
axes.set_aspect(1)
# plt.scatter(pointings_vec[:,0], pointings_vec[:,1], s=1, c='r', label='Pointings')
plt.scatter(plot_targets.ra.degree, plot_targets.dec.degree, s=0.1, c='b', label='Targets', alpha = 1)

for i in range(0, len(pointings_vec)):
    circle120 = plt.Circle((pointings_vec[i]), 1.295, fill = False, lw = 0.15, color = 'green')
    axes.add_artist(circle120)

plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
plt.savefig('plots/v3_gaia_query.png', dpi=300)

# --- Beam Coordinates --- 
beam_ra = []; beam_dec = []
for i in range(0, len(sources)): 
    beam_ra.append(pointings_ra[in_beam_targets[1][i]])
    beam_dec.append(pointings_dec[in_beam_targets[1][i]])

# --- Saving ids to file ---
np.savetxt('txt/gaia_ids.txt', np.unique(sources), fmt='%s')
gaia_df = pd.DataFrame({'source_id': sources, 'ra': plot_targets.ra.degree, 'dec': plot_targets.dec.degree, 'seperation': in_beam_targets[2].degree, 'beam': in_beam_targets[1] , 'beam_ra': beam_ra, 'beam_dec': beam_dec})
gaia_df.to_csv('csv/gaia_query.csv', index=False)

script_runtime_end = time.time()
print('Time taken for script to run %s secs.' % ("{:.1f}".format(script_runtime_end - script_runtime_start)))