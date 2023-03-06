""" 
CODE FUNCTION: To plot .h5 files into dynamic waterfall plots using I-LOFAR and LOFAR-SE data that is a product of the Breakthrough Listen pipelines. 
Authors: Charles Giese, Owen Johnson
Date last updated: 22/07/2022
"""
from blimpy import Waterfall
import numpy as np
import matplotlib.pyplot as plt
import os
from turbo_seti.find_event.find_event import read_dat
# from turbo_seti.find_event.plot_event import plot_waterfall
from astropy.time import Time
import matplotlib
from blimpy.utils import rebin
import sys

fontsize=16
font = {'family' : 'serif', 'serif':'cmr10', 'size' : fontsize} 
# font = {'family' : 'DejaVu Sans',
# 'size' : fontsize}
MAX_IMSHOW_POINTS = (4096, 1268)


#makes plot of every hit in a .dat file

try:
	os.mkdir('./hit_pngs')
except:
        print('Could not create directory')


def overlay_drift(f_event, f_start, f_stop, drift_rate, t_duration, offset='auto'):
    r'''
    Creates a dashed red line at the recorded frequency and drift rate of
    the plotted event - can overlay the signal exactly or be offset by
    some amount (offset can be 0 or 'auto').
    '''
    # determines automatic offset and plots offset lines
    if offset == 'auto':
        offset = ((f_start - f_stop) / 10)
        plt.plot((f_event - offset, f_event),
                 (10, 10),
                 "o-",
                 c='#cc0000',
                 lw=2)

    # plots drift overlay line, with offset if desired
    plt.plot((f_event + offset, f_event + drift_rate/1e6 * t_duration + offset),
             (0, t_duration),
             c='#cc0000',
             ls='dashed', lw=2)

filename=sys.argv[1]+'.h5'
dat_file=sys.argv[1]+'.dat'
on_source_name= filename.split('.')[0]

df=read_dat(dat_file)
low_freqs=df['FreqStart']
high_freqs=df['FreqEnd']
snr=df['SNR']
drift_rates=df['DriftRate']

matplotlib.rc('font', **font)

for i in range(len(low_freqs)):
	f_start=low_freqs[i]
	f_stop=high_freqs[i]

	plt.figure(i)
	fil = Waterfall(filename, f_start=f_start, f_stop=f_stop)

	dummy, plot_data = fil.grab_data()

# rebin data to plot correctly with fewer points
	dec_fac_x, dec_fac_y = 1, 1
	if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
		dec_fac_x = plot_data.shape[0] / MAX_IMSHOW_POINTS[0]
	if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
		dec_fac_y =  int(np.ceil(plot_data.shape[1] /  MAX_IMSHOW_POINTS[1]))
	plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)


	f_mid = round(np.abs(f_start+f_stop)/2.,  4)
	mid_f = round(np.abs(f_start+f_stop)/2.,  4)
	drift_rate=drift_rates[i]

	# read in data

	fig=plt.figure(1)

	t0 = fil.header['tstart']
	# make plot with plot_waterfall
	source_name = 'B2217+47'

	# calculate parameters for estimated drift line
	t_elapsed = Time(fil.header['tstart'], format='mjd').unix - Time(t0, format='mjd').unix
	t_duration = (fil.n_ints_in_file - 1) * fil.header['tsamp']
	f_event = f_mid + drift_rate / 1e6 * t_elapsed
	offset= 0

	# calculate the width of the plot based on making sure the full drift is visible
	bandwidth = 2.4 * abs(drift_rate)/1e6 * t_duration
	bandwidth = np.max((bandwidth, 500./1e6))

	# Get start and stop frequencies based on midpoint and bandwidth
	f_start, f_stop = np.sort((f_mid - (bandwidth/2),  f_mid + (bandwidth/2)))

	this_plot = fil.plot_waterfall(f_start=f_start, f_stop=f_stop)

    # plot estimated drift line
	overlay_drift(f_event, f_start, f_stop, drift_rate, t_duration, offset)

	# Title the full plot

	plot_title = "%s \n $\dot{\\nu}$ = %2.3f Hz/s , MJD:%5.5f" % (on_source_name, drift_rate, t0)

	plt.title(plot_title)
	# Format full plot
	plt.xticks(np.linspace(f_start, f_stop, num=4), ['','','',''])

	# More overall plot formatting, axis labelling
	factor = 1e6
	units = 'Hz'

	ax = plt.gca()

	xloc = np.linspace(f_start, f_stop, 5)
	xticks = [round(loc_freq) for loc_freq in (xloc - mid_f)*factor]
	if np.max(xticks) > 1000:
		xticks = [xt/1000 for xt in xticks]
		units = 'Hz'
	plt.xticks(xloc, xticks)
	plt.xlabel("Relative Frequency [%s] from %f MHz"%(units,mid_f),fontdict=font)

	plt.savefig('./hit_pngs/'+'hit_num'+str(i)+'.png', bbox_inches='tight')
	plt.close('all')


