"""
Code Purpose: To process barycentrically correct 0000.fil files in a specified directory.
Author: Owen Johnson
Last Major Update: 25/10/2022
Notes: Calls on code created by Dr. Vishal Gajjar, see https://github.com/gajjarv/BaryCentricCorrection for more information
""" 

import argparse
import os 
import re
import glob
import pandas as pd 
import numpy as np
from datetime import datetime
from email_operator import email_users
from astropy import units as u
from astropy.coordinates import SkyCoord

'''
--- FUNCTIONS --- 
'''

def coord_bary(coords, coordinate_idx, index): # Puts ra and dec into format for barycentering. Co-ordinate in astropy, index [0 or 1] for RA or DEC respectively. 
    bcord_idx = coords.to_string('hmsdms')[coordinate_idx].split(' ')[index] 
    bcord_flt = re.sub("[^0-9]", "", str(bcord_idx))
    bcord_new = bcord_flt[:6] + '.' + bcord_flt[6:]
    return bcord_new



# - Command line arguments - 
parser = argparse.ArgumentParser(description="Processes LOFAR filterbanks through Vishal Gajjar's Barycentering Code.")
parser.add_argument('-i', '--path',  help="path to input .fil files.")
parser.add_argument('-s', '--station', help="Which station is this data being processed on? Irl or Swe")
args = parser.parse_args()

# - setting station specific paths - 
station_str = str(args.station.lower()).split(' ')[0]

if station_str == 'sweden':
    bary_path = '/home/obs/linux_64/BaryCentricCorrection/barycentre_seti'
    station_tag = 'p'
if station_str == 'ireland':
    bary_path = '/home/owen/BaryCentricCorrection/barycentre_seti'
    station_tag = 'n'
else: 
    raise ValueError('Please select a valid station.')

# - Loading file paths into an array - 
input_data_directory = args.path
fil_list = sorted(glob.glob(input_data_directory + '/TIC*0000.fil'))

# - To find the file location no matter how deep the path is - 
sample_path = fil_list[0]
if 'TIC' in sample_path:
    idx = sample_path.index('TIC')

# - Extracting the TIC id from the path name - 
target_ids = []
for file_path in fil_list:
    target_ids.append(int(file_path[idx + 3:].split('.')[0]))
target_ids = np.array(target_ids)

# - Finding target parameters from the masterlist of targets of interest - 
df_mask = []
ext_pd = pd.read_csv('TESS_ext_target_data.csv')
for tic_id in target_ids: 
     df_mask.append(int(ext_pd[ext_pd['TIC_ID']==tic_id].index.values))

ext_pd_masked = ext_pd.loc[df_mask] # - masking the dataframe. 

coords = SkyCoord(ra=ext_pd_masked['ra']*u.degree, dec=ext_pd_masked['dec']*u.degree, frame='icrs')

for idx in range(0, len(coords)): 
    ra = coord_bary(coords, idx, 0)
    dec = coord_bary(coords, idx, 1)

    fil_output  = '%s/TIC%s.bary.0000.fil' % (args.path, str(target_ids[idx]))

    os.system('filedit --ra %s --dec %s --src-name TIC-Obj %s' % (ra, dec, fil_list[idx])) # - Ensuring header has the correct ra and dec values 
    os.system("%s %s -site %s -verbose > %s" % (bary_path, fil_list[idx], station_tag, fil_output)) # - Calling the barycenter code 
