"""
Code Purpose: Slice filterbanks every 10 Mhz to make it more digestible to use in Breakthrough Listen analysis pipelines. 
Author: Owen Johnson.
Last Major Update: 03/08/2022
"""

import numpy as np 
import argparse
import glob 
import os 
from tqdm import tqdm # Used for the loading bar seen in the terminal. 
# from multiprocessing import Pool # TO DO: Add multi-processing 

# - Command line arguments - 
parser = argparse.ArgumentParser(description='Slices the filterbanks up for processing.')
parser.add_argument('-i', '--path',  help="path to input .fil files.")
parser.add_argument('-o', '--output',  help="path to output for diced .fil files.")
args = parser.parse_args()

# - Loads all .fil files in a certain path into an array of strings - 
file_list = sorted(glob.glob(str(args.path) + '*.fil'))

intervals = 9 # - Hard coded so the files are split every 10 Mhz, easily changed for desired frequency width. 
freq_rngs = np.linspace(110, 190, intervals)

def main(file_idx): 
    for idx in range(0, (intervals - 1)): 
        output_name = file_list[file_idx].split('/')[-1]
        output_path = args.output + 'sliced.' +  str(round(freq_rngs[idx])) + '.' + str(round(freq_rngs[idx + 1])) + '.' + output_name 
        os.system("bldice -f %s -o %s -b %s -e %s" % (file_list[file_idx], output_path, freq_rngs[idx], freq_rngs[idx + 1])) # - calling bldice using os as it's a command line util. 

for i in tqdm(range(0, len(file_list))): # - Slicing for each file in the path directory. 
    main(i)

print('Slicing Done!')