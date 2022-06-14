#! /usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

parser = argparse.ArgumentParser(description='Compares .dat files from turboSETI searches ')
parser.add_argument("-i", "--inputfiles", nargs='+', help="paths to input .dat files. Irish station comes first then Swedish.")
parser.add_argument("-s", "--snr_cutoff", help = "Maximum SNR threshold")
args = parser.parse_args()
ie_dat = args.inputfiles[0]
se_dat = args.inputfiles[1]
snr_thresh = args.snr_cutoff

print('Irish file path:', ie_dat, 'Shape:', np.shape(ie_dat))
print('Swedish file path:', se_dat, 'Shape:', np.shape(se_dat))
print('Inputted SNR threshold:', np.array(snr_thresh), 'Type:', type(snr_thresh))

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

df_ie = read_dat('./'+ie_dat)
df_se = read_dat('./'+se_dat)

dr = 0.1 #Hz/s
df = 5 #Hz

hit_matches = []

for ind in df_ie.index:
	hit_num = df_ie['TopHitNum'][ind]
	frequency = df_ie['Freq'][ind] * 1e6
	drift_rate = df_ie['DriftRate'][ind]
	sigma = df_ie['SNR'][ind]

	for ind_s in df_se.index:
		hit_num_se = df_se['TopHitNum'][ind_s]
		frequency_se = df_se['Freq'][ind_s] *1e6
		drift_rate_se = df_se['DriftRate'][ind_s]
		sigma_se = df_se['SNR'][ind_s]
		if sigma < snr_thresh:
			if sigma_se < snr_thresh:
				if frequency - df <= frequency_se <= frequency + df:
					if drift_rate - dr <= drift_rate_se <= drift_rate +dr:
						hit_matches.append('IE'+str(hit_num)+'_SE'+str(hit_num_se))
						print('Hit Numbers: IR ', ind, 'SE ', ind_s)
						print('IE Freq: ',frequency, 'SE Freq: ', frequency_se)
						print('IE D_rate:', drift_rate, 'SE D_rate: ', drift_rate_se)
