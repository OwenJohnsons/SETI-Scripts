import glob 
import argparse
import numpy as np
import plotille

# --- CLI Arguments ---
parser = argparse.ArgumentParser(description='Long term statistics')
parser.add_argument('-p', '--path', type=str, required=True, help='Path to data')

path = parser.parse_args().path

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

# --- Fetching .dat files --- 
dat_files = glob.glob(path + '/**/*.dat', recursive=True)
fil_files = glob.glob(path + '/**/TIC*bary*.fil', recursive=True)

print('==========================================')
print('Number of .dat files found:', len(dat_files))
print('Number of barycentered .fil files found:', len(fil_files))
print('==========================================')
# seperating tic ids
ids = []
for file in dat_files:
    tic_id = file.split('/')[-1].split('.')[0]
    ids.append(tic_id)
unique, unique_indices = np.unique(ids, return_index=True)

if len(unique_indices) != len(ids):
    duplicates = np.delete(ids, unique_indices)
    print('WARNING: %s duplicate TIC IDs found.' % len(duplicates))
    
    
drift_rates = []; snrs = []; corrected_freqs = []
for file in dat_files: 
    data = dat_file_reader(file)
    drift_rates = np.append(drift_rates, data[0])
    snrs = np.append(snrs, data[1])
    corrected_freqs = np.append(corrected_freqs, data[2])
    
print('Number of hits detected across %s targets:' % len(unique_indices), len(drift_rates))
print('---')
print('Drift Rate Statistics:')
print('Mean:', np.mean(drift_rates))
print('Median:', np.median(drift_rates))
print('Standard Deviation:', np.std(drift_rates))
print('---')
print('SNR Statistics:')
print('Mean:', np.mean(snrs))
print('Median:', np.median(snrs))
print('Standard Deviation:', np.std(snrs))
print('---')
print('Corrected Frequency Statistics:')
print('Mean:', np.mean(corrected_freqs))
print('Median:', np.median(corrected_freqs))
print('Standard Deviation:', np.std(corrected_freqs))
print('---')
print(plotille.histogram(
    drift_rates,
    bins=100,
    width=60,
    height=30,
    X_label='Drift Rate [Hz/s]',
    Y_label='Counts',
    linesep='\n'
))
