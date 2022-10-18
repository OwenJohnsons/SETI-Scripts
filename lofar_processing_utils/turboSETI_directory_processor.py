"""
Code Purpose: To process filterbanks in an active directory
Author: Owen Johnson
Last Major Update: 14/10/2022
Notes: Uses pieces of code from BL internship matierals. 
""" 

import argparse
import os 
import glob
from tqdm import tqdm 
from turbo_seti.find_doppler.find_doppler import FindDoppler
from datetime import datetime
from email_operator.py import email_users
import numpy as np 


# - Command line arguments - 
parser = argparse.ArgumentParser(description='Processes LOFAR filterbanks through turboSETI.')
parser.add_argument('-i', '--path',  help="path to input .fil files.")
parser.add_argument('-t', '--processed_trgts', help = "A .txt file with an array of previously processed targets such that the same targets aren't processed twice")
args = parser.parse_args()

input_data_directory = args.path 
processed_trgts = np.loadtxt(args.processed_trgts, dtype='str')

try:
	os.mkdir(input_data_directory + '/dat_files')
except:
        print('Could not create directory')

output_data_directory = input_data_directory + '/dat_files'
fil_list = glob.glob(input_data_directory + '/TIC*.fil')

for file in tqdm(fil_list):
    if file not in processed_trgts:
        try: 
            doppler = FindDoppler(file,
                            max_drift = 4, # Max drift rate = 4 Hz/second
                            snr = 10,      # Minimum signal to noise ratio = 10:1
                            n_coarse_chan = 411,
                            out_dir = output_data_directory, # This is where the turboSETI output files will be stored.
                            gpu_id = 0 # Uses default gpu on each LOFAR station. 
                            )
            doppler.search()
            fil_list.append(file)
            np.savetxt(args.processed_trgts, fil_list,  fmt="%s")

        except: 
            print('Turbo search file failed on this target... moving to next file')
            failed_trgts = open(output_data_directory + "/failed_trgts.txt","w+")
            time_stamp = datetime.now()
            failed_trgts.write("/n", time_stamp, file)
            failed_trgts.close()
    else: 
        print('Critical Error, check .log')
    
print('\n Job Finished')
email_users()