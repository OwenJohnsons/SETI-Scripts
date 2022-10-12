"""
Code Purpose: Quickly Plot .h5 files when testing SETI-gen on the LOFAR pipeline. 
Author: Owen Johnson.
Last Major Update: 04/08/2022
"""

import blimpy as bl
from blimpy import Waterfall
import glob
import matplotlib.pyplot as plt

plt.rc('font', family = 'serif', serif = 'cmr10') 
plt.rcParams.update({'font.size': 22})
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'stix'

files = glob.glob('*.h5')

for file in files:
    h5 = Waterfall(file)
    plt.figure(figsize=(10,8))
    h5.plot_waterfall()
    plt.savefig(str(file.split('.')[0] + '.png'),  bbox_inches='tight')