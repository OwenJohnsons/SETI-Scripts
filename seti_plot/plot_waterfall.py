''' 
Code Purpose: Plots a inputted waterfall file at a specified central frequency. 
Author: Owen A. Johnson 
Date of last major modification: 06/03/2023
'''
import blimpy as bl 
from blimpy import Waterfall
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import argparse 
import os 

# --- Command line arguments ---
parser = argparse.ArgumentParser(description='Plots a waterfall file given and input and a central frequency.')
parser.add_argument('--filename', '-i', type=str, help='Path of the file to be plotted.')
parser.add_argument('--central_freq', '-f', type=float, help='Central frequency (MHz) of the file to be plotted.')
args = parser.parse_args()

# - loading filterbank file and details -
filename = args.filename
central_frq = args.central_freq
start_frq = (central_frq) - 250/1e6; end_frq = (central_frq) + 250/1e6

source_name = filename.split('/')[-1].split('.')[0]
path = filename[0:-len(filename.split('/')[-1])]

try :
    os.mkdir('%sblimpy_plots' % path)
except :
    print('Could not create directory.')
    
save_path = '%sblimpy_plots/%s' % (path, (source_name + '_' + str(central_frq)))

# --- Loading Filterbank ---
print('Bandwidth to be plotted ', np.round((end_frq - start_frq)*1e6, 5), ' Hz')
data = Waterfall(args.filename, load_data=True, f_start = start_frq, f_stop = end_frq)

# --- Plot Parameters --- 
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2

# --- Plotting Filterbank --- 
fig = plt.figure(figsize=(10, 10))
data.plot_spectrum(logged=True,  f_start=start_frq, f_stop=end_frq)
plt.title(source_name)
plt.savefig('%s_spectrum.png' % save_path, dpi=300)

fig = plt.figure(figsize=(14, 10))
data.plot_waterfall(f_start=start_frq, f_stop=end_frq, cmap = 'cool')
plt.title(source_name)

xloc = np.linspace(start_frq, end_frq, 5)
xticks = [round(loc_freq) for loc_freq in (xloc - central_frq)*1e6]; units = 'Hz'
plt.xticks(xloc, xticks)
plt.xlabel("Relative Frequency [%s] from %s MHz"%(units, central_frq))
plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
plt.savefig('%s_waterfall.png' % save_path, dpi=300, bbox_inches='tight')