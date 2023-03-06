''' 
Code Purpose: Plots a inputted waterfall file. 
Author: Owen A. Johnson 
Date of last major modification: 
'''

from blimpy import Waterfall
import numpy as np
import matplotlib.pyplot as plt
import argparse 
import os 

# --- Command line arguments ---
parser = argparse.ArgumentParser(description='Plots a waterfall file given and input and a central frequency.')
parser.add_argument('--filename', '-i', type=str, help='Path of the file to be plotted.')
parser.add_argument('--central_freq', '-f', type=float, help='Central frequency (MHz) of the file to be plotted.')
args = parser.parse_args()

filename = args.filename
central_frq = args.central_freq
start_frq = (central_frq) - 250/1e6; end_frq = (central_frq)*1e6 + 250/1e6

source_name = filename.split('/')[-1].split('.')[0]
path = filename[0:-len(filename.split('/')[-1])]

try :
    os.mkdir('%sblimpy_plots' % path)
except :
    print('Could not create directory.')
    
save_path = '%sblimpy_plots/%s' % (path, source_name)

# --- Plotting Filterbank --- 
fig = plt.figure(figsize=(10, 10))
obs = Waterfall(filename, f_start=start_frq, f_stop=end_frq)

fig = plt.figure(figsize=(10, 10))
obs.plot_spectrum(logged=True,  f_start=start_frq, f_stop=end_frq)
plt.title(source_name)
plt.savefig('%s_spectrum.png' % save_path, dpi=300)

fig = plt.figure(figsize=(10, 10))

dummy, plot_data = obs.grab_data()

# rebin data to plot correctly with fewer points
dec_fac_x, dec_fac_y = 1, 1
if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
    dec_fac_x = plot_data.shape[0] / MAX_IMSHOW_POINTS[0]
if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
    dec_fac_y =  int(np.ceil(plot_data.shape[1] /  MAX_IMSHOW_POINTS[1]))
plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)

obs.plot_waterfall(cmap='viridis',f_start=start_frq, f_stop=end_frq)
plt.title(source_name)
plt.xticks(np.linspace(start_frq, end_frq, num=4), ['','','',''])
plt.savefig('%s_waterfall.png' % save_path, dpi=300)