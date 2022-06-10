#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib
from sigpyproc.Readers import readTim
from astropy.time import Time

parser = argparse.ArgumentParser(description='Plot .tim timeseries files')
parser.add_argument("-i", "--inputfiles", nargs='+', help=".tim files, IR then SE")
parser.add_argument("-s", "--time_start", help="Start time for plotting")
parser.add_argument("-e", "--time_end", help="End time for plotting")
parser.add_argument("-o", "--outputfile", help="Name of output .png")
args = parser.parse_args()
inputfiles = args.inputfiles
time_start = float(args.time_start)
time_end = float(args.time_end)
outputfile = args.outputfile

fontsize=12
font = {'family' : 'DejaVu Sans',
'size' : fontsize}
matplotlib.rc('font', **font)

tim_IR = readTim(inputfiles[0])
tim_SE = readTim(inputfiles[1])

#trying to figure out time stuff:

source_name = tim_IR.header.source_name

obs_duration_IR = tim_IR.header.nsamples * tim_IR.header.tsamp
tim_flux_IR = np.array(tim_IR)
time_IR = np.linspace(0, obs_duration_IR, np.size(tim_flux_IR))

obs_duration_SE = tim_SE.header.nsamples * tim_SE.header.tsamp
temp = np.array(tim_SE)
tim_flux_SE = temp[int(17.86/tim_SE.header.tsamp):]
time_SE = np.linspace(0, obs_duration_SE, np.size(tim_flux_SE))



fig = plt.figure(1)
ax_IR = fig.add_subplot(211)
ax_IR.plot(time_IR, tim_flux_IR, c='k')
ax_IR.set_xlim(time_start, time_end)
ax_IR.set_xlabel('Time (s)')
ax_IR.set_title(source_name+' IR')
ax_IR.set_ylabel('Flux')
ax_SE = fig.add_subplot(212)
ax_SE.plot(time_SE, tim_flux_SE, c='k')
ax_SE.set_xlim(time_start, time_end)
ax_SE.set_xlabel('Time (s)')
ax_SE.set_title(source_name+' SE')
ax_SE.set_ylabel('Flux')
plt.tight_layout()
plt.savefig(source_name+'_'+str(time_start)+'-'+str(time_end)+'sec.png')
plt.show()
