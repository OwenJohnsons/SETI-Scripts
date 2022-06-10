from blimpy import Waterfall
from blimpy import dice
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os
import argparse
from turbo_seti.find_event.find_event import read_dat
from turbo_seti import plot_event

parser = argparse.ArgumentParser(description='Plots hits in a very rough and ready manner, not super useful right now')

parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Name of .fil/.dat file")
parser.add_argument("-o", "--out_dir", help="Directory to save images to")
#parser.add_argument("-fl", "--lower_fbound", help="lower bound of frequency range")
#parser.add_argument("-fu", "--upper_fbound", help="upper bound of frequency range")


args = parser.parse_args()
inputfile=args.inputfile
out_dir=args.out_dir

fil_file=inputfile+'.h5'
dat_file=inputfile+'.dat'

df=read_dat(dat_file)
low_freqs=df['FreqStart']
high_freqs=df['FreqEnd']
snr=df['SNR']
drift_rates=df['DriftRate']



for i in range(len(low_freqs)):
	f0=low_freqs[i]
	f1=high_freqs[i]


	print('Loading data in range:', f0, f1)
	obs=Waterfall(fil_file, f_start=f0, f_stop=f1)

	fig=plt.figure(i)
	obs.plot_waterfall()
	plt.savefig(out_dir+str(f0)+'_'+str(f1)+'_SNR'+str(snr[i])+'.png')
	fig.close()
