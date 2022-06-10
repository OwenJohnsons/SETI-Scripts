#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib
from sigpyproc.Readers import readTim


parser = argparse.ArgumentParser(description='Timeseries plotting')
parser.add_argument("-i", "--inputfile", help=".time file")
parser.add_argument("-s", "--time_start", help="Start point for plotting")
parser.add_argument("-e", "--time_end", help="End point for plotting")
args = parser.parse_args()
inputfile = args.inputfile
time_start = args.time_start
time_end = args.time_end
outputfile = args.outputfile

fontsize=12
font = {'family' : 'DejaVu Sans',
'size' : fontsize}
matplotlib.rc('font', **font)

tim = readTim(inputfile)

source_name = tim.header.source_name
obs_duration = tim.header.nsamples * tim.header.tsamp

tim_flux = np.array(tim)
time = np.linspace(0, obs_duration, np.size(tim_flux))


fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(time, tim_flux, c='k')
ax.set_xlim(time_start, time_end)
ax.set_xlabel('Time (s)')
ax.set_title(source_name)
ax.set_ylabel('Flux')
plt.show()
