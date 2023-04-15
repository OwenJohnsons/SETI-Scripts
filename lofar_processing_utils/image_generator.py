'''
Code Purpose: Takes filterbanks and turboSETI search outputs and generates plots of the hits with a 500 Hz frequency range around the corrected frequency. This is an optional part of the overall pipeline. An average filterbank with ~300 hits across the 90 MHz band takes X hours to generate.
Author: Owen A. Johnson 
Date of last edit: 15/04/2023
'''

import argparse
import os, sys
import numpy as np
import glob 
import matplotlib.pyplot as plt
import smplotlib 
import blimpy as bl
from tqdm import tqdm


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

# --- Command line arguments ---
parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', type=str, help='Path to directory containing filterbank and .dat files.', required=True)
parser.add_argument('-v', '--verbose', type=bool, help='Prints out extra information for debugging.', required=False, default=False)
parser.add_argument('-s', '--station', type=str, help='LOFAR Station: IRE or SWE.', required=True)
parser.add_argument('-spectra', '--spectra', type=bool, help='Generate spectra plots.', required=False, default=False)

path = parser.parse_args().path
verbose = parser.parse_args().verbose
station = parser.parse_args().station

if station not in ['IRE', 'SWE']:
    raise ValueError('Please enter a valid station name i.e. IRE or SWE.')

if station == 'IRE':
    station = 'I-LOFAR'
elif station == 'SWE':
    station = 'LOFAR-SE'

# --- Finding File Paths ---
filterbanks = glob.glob(path + '/**/*.fil', recursive=True)
dat_files = glob.glob(path + '/**/*.dat', recursive=True)

'''
.dat file headers are as follows:
1 - Drift Rate 
2 - SNR 
3 - Uncorrected Frequency [MHz]
4 - Corrected Frequency [MHz]
6 - Start Frequency [MHz]
7 - End Frequency [MHz]
'''
def dat_file_reader(file):
    contents = np.loadtxt(file, skiprows = 9, usecols = (1, 2, 4, 6, 7))
    return contents.T

os.makedirs(path + 'plots', exist_ok=True)

# --- See if filterbank files have a matching .dat file ---
with tqdm(total=len(filterbanks), position=0, leave=True) as pbar:
    for file in tqdm(filterbanks, position=0, leave=True):
        tic_id = file.split('/')[-1].split('.')[0][3:]
        for dat_file in dat_files:
            if tic_id in dat_file:
                if verbose:
                    print('Match Found')
                    print(file)
                    print(dat_file)
                    print('')

                # --- Making output directory for plots ---
                output_path = path + 'plots/TIC' + tic_id
                os.makedirs(output_path, exist_ok=True)

                hit_details = dat_file_reader(dat_file)
                width = 250/1e6
                
                for i in range(0, len(hit_details[0])):
                    blockPrint()

                    corrected_freq = hit_details[2][i]; drift_rate = hit_details[0][i]
                    filterbank = bl.Waterfall(file, f_start = corrected_freq - width, f_stop = corrected_freq + width) # - Taking a frquency subset 

                    # --- Tick Setup ---
                    xloc = np.linspace(corrected_freq - width, corrected_freq + width, 5)
                    xticks = [round(loc_freq) for loc_freq in (xloc - corrected_freq)*1e6]; units = 'Hz'

                    # === Plotting Spectra === 
                    if parser.parse_args().spectra:
                        plt.figure(figsize=(5, 5))
                        filterbank.plot_spectrum()
                        plt.title('TIC%s (%s)\nDrift Rate: %s Hz/s' % (tic_id, station, np.round(hit_details[0][i], 3)))
                        plt.xticks(xloc, xticks)
                        plt.xlabel("Relative Frequency [%s] from %s MHz"%(units, np.round(corrected_freq, 3)))
                        plt.savefig(output_path + '/TIC%s_hit_%s_spectrum.png' % (tic_id, i + 1), dpi = 200)
                        plt.tight_layout()
                        plt.show()

                    # === Plotting Waterfall ===
                    plt.figure(figsize=(5, 6))
                    filterbank.plot_waterfall()
                    plt.title('TIC%s (%s)\nDrift Rate: %s Hz/s' % (tic_id, station, np.round(hit_details[0][i], 3)))
                    plt.xticks(xloc, xticks)
                    plt.xlabel("Relative Frequency [%s] from %s MHz"%(units, np.round(corrected_freq, 3)))
                    plt.tight_layout()
                    plt.savefig(output_path + '/TIC%s_hit_%s_waterfall.png' % (tic_id, i + 1), dpi = 200)
                    plt.show()

                    enablePrint()