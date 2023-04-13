'''
Code Purpose: Given a path will find all filterbank files in a directory, change the ra and dec in their header and then convert them using Vishal Gajjar's barycentering code before processing.  
Author: Owen A. Johnson 
Date of last major change: 12/04/2023
'''

import argparse
import os 
import glob
import pandas as pd 
import astropy.units as u
from astropy.coordinates import SkyCoord

# --- Command line arguments ---
parser = argparse.ArgumentParser(description='Barycentering automation')
parser.add_argument('-p', '--path', type=str, help='Path to directory containing filterbank files.', required=True)
parser.add_argument('-o', '--output', type=str, help='Path to directory where barycentered files will be stored.', required=True)
parser.add_argument('-s', '--station', type=str, help='''Station where data is being processed. Key : Station | IRE : I-LOFAR | SWE : LOFAR-SE |''', required=True)

path = parser.parse_args().path
output = parser.parse_args().output
station = parser.parse_args().station

# --- Setting Station --- 
if station == 'SWE':
    station = 'p'
    software_path = '/home/obs/linux_64/BaryCentricCorrection/barycentre_seti'
    print('Station set to LOFAR-SE')
elif station == 'IRE':
    station = 'n'
    software_path = '/home/owen/BaryCentricCorrection/barycentre_seti'
    print('Station set to I-LOFAR')

print('Number of files to be barycentered: %s' % len(glob.glob(path + '/**/TIC*rawspec.0000.fil', recursive=True)))

# --- Creating Output Directory if none ---
if not os.path.exists(output):
    os.makedirs(output)
    
# --- Finding File Paths ---
files = glob.glob(path + '/**/TIC*rawspec.0000.fil', recursive=True)
TESS_csv = pd.read_csv('/home/owen/SETI-Scripts/observation_lists/TESS_ext_target_data.csv')

for file in files: 
    sub_dir = file.split('/')[-1]
    tic_id = sub_dir.split('.')[0][3:]
    trgt_data = TESS_csv[TESS_csv['TIC_ID'] == int(tic_id)] # - isolating target data from file title. 
    # - Putting Coordinates in correct format for barycentering code. - 
    ra = trgt_data['ra'].values[0]; dec = trgt_data['dec'].values[0]
    c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    ra = c.ra.to_string(unit=u.hour, sep='', pad=True, precision=3); dec = c.dec.to_string(sep='', pad=True, alwayssign=True, precision=3)
    # - Changing header values -
    os.system('filedit --ra %s --dec %s --src-name TIC-Obj %s' % (ra, dec, file))

    # - Barycentering -
    o_filename = file.split('/')[-1].replace('rawspec', 'bary')
    print(o_filename)
    os.system('%s %s -site %s -verbose > %s' % (software_path, file, station, (output + '/' + o_filename)))


