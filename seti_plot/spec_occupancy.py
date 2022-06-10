#! /usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys
import time
import re
import pandas as pd

#script for plotting spectral and drift rate occupancy of hits

fontsize=16
font = {'family' : 'DejaVu Sans',
'size' : fontsize}
MAX_IMSHOW_POINTS = (4096, 1268)


def read_dat(filename):
    r"""
    Read a turboseti .dat file.
    Parameters
    ----------
    filename : str
        Name of .dat file to open.
    Returns
    -------
    df_data : dict
        Pandas dataframe of hits.
    """
    file_dat = open(filename.strip())
    hits = file_dat.readlines()

    # Get info from the .dat file header
    FileID = hits[1].strip().split(':')[-1].strip()
    Source = hits[3].strip().split(':')[-1].strip()

    MJD = hits[4].strip().split('\t')[0].split(':')[-1].strip()
    RA = hits[4].strip().split('\t')[1].split(':')[-1].strip()
    DEC = hits[4].strip().split('\t')[2].split(':')[-1].strip()

    DELTAT = hits[5].strip().split('\t')[0].split(':')[-1].strip()  # s
    DELTAF = hits[5].strip().split('\t')[1].split(':')[-1].strip()  # Hz

    # Get info from individual hits (the body of the .dat file)
    all_hits = []
    for hit_line in hits[9:]:
        hit_fields = re.split(r'\s+', re.sub(r'[\t]', ' ', hit_line).strip())
        all_hits.append(hit_fields)

    # Now reorganize that info to be grouped by column (parameter)
    # not row (individual hit)
    if all_hits:
        TopHitNum = list(zip(*all_hits))[0]
        DriftRate = [float(df) for df in list(zip(*all_hits))[1]]
        SNR = [float(ss) for ss in list(zip(*all_hits))[2]]
        Freq = [float(ff) for ff in list(zip(*all_hits))[3]]
        ChanIndx = list(zip(*all_hits))[5]
        FreqStart = list(zip(*all_hits))[6]
        FreqEnd = list(zip(*all_hits))[7]
        CoarseChanNum = list(zip(*all_hits))[10]
        FullNumHitsInRange = list(zip(*all_hits))[11]

        data = {'TopHitNum': TopHitNum,
                'DriftRate': DriftRate,
                'SNR': SNR,
                'Freq': Freq,
                'ChanIndx': ChanIndx,
                'FreqStart': FreqStart,
                'FreqEnd': FreqEnd,
                'CoarseChanNum': CoarseChanNum,
                'FullNumHitsInRange': FullNumHitsInRange
                }

        # Creating pandas dataframe from data we just read in
        df_data = pd.DataFrame(data)
        df_data = df_data.apply(pd.to_numeric)

    else:
        df_data = pd.DataFrame()

    # Matching column information from before to the .dat data we read in
    df_data['FileID'] = FileID
    df_data['Source'] = Source.upper()
    df_data['MJD'] = MJD
    df_data['RA'] = RA
    df_data['DEC'] = DEC
    df_data['DELTAT'] = DELTAT
    df_data['DELTAF'] = DELTAF

    # Adding extra columns that will be filled out by this program
    df_data['Hit_ID'] = ''
    df_data['status'] = ''
    df_data['in_n_ons'] = ''
    df_data['RFI_in_range'] = ''

    return df_data

dat_file = sys.argv[1]
on_source_name = 'B2217+47'

df = read_dat(dat_file)

freqs = df['Freq'].values
bins = np.linspace(np.floor(np.min(freqs)), np.floor(np.max(freqs)), int(np.floor(np.max(freqs)) - np.floor(np.min(freqs)) +1))

drift_rates = df['DriftRate'].values
drift_bins = 10000

print(drift_rates)

fig = plt.figure(1)
ax = fig.add_subplot(121)
ax.hist(freqs, np.size(bins), edgecolor='k')
ax.set_xlabel('Frequency (MHz)', fontdict=font)
ax.set_ylabel('Hit count', fontdict=font)

ax2 = fig.add_subplot(122)
ax2.hist(drift_rates, 100, edgecolor='k')
ax2.set_xlabel('Drift rate, Hz/s', fontdict=font)
ax2.set_ylabel('Hit count', fontdict=font)
ax2.set_yscale('log')
plt.tight_layout()
plt.show()
