"""
Code Purpose: To process barycentrically correct 0000.fil files in a specified directory.
Author: Owen Johnson
Last Major Update: 18/10/2022
Notes: Calls on code created by Dr. Vishal Gajjar, see https://github.com/gajjarv/BaryCentricCorrection for more information
""" 

import argparse
import os 
import glob
import pandas as pd 
import numpy as np
from datetime import datetime
from email_operator import email_users


  
# - Command line arguments - 
parser = argparse.ArgumentParser(description="Processes LOFAR filterbanks through Vishal Gajjar's Barycentering Code.")
parser.add_argument('-i', '--path',  help="path to input .fil files.")
args = parser.parse_args()

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