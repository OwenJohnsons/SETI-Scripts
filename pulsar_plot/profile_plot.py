#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib

parser = argparse.ArgumentParser(description='Plot pulse profile')
parser.add_argument("-i", "--inputfile", help="Pulse profile file")
parser.add_argument("-s", "--source_name", help="Name of Source")
parser.add_argument("-o", "--outputfile", help="Name of output")
args = parser.parse_args()
inputfile = args.inputfile
source = args.source_name
outputfile = args.outputfile

fontsize=12
font = {'family' : 'DejaVu Sans',
'size' : fontsize}
matplotlib.rc('font', **font)

names = ["bin", "amp"]
data = pd.read_csv(inputfile, sep='\s+', skiprows=1, names=names)
bin = data["bin"].values
amp = data["amp"].values

bins = bin / np.max(bin)
flux = -(amp / np.min(amp))

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(bins, flux, c='k')
ax.set_xlabel('Pulse Phase')
ax.set_ylabel('Flux')
ax.set_title(source+' Folded Pulse Profile')
plt.savefig(outputfile+'.png')
