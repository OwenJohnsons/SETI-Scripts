'''
Code Purpose: Makes a table of extended TESS Target info from TESS inputed TESS IDs 
Author: Owen Johnson 
Last Major Update: 17/01/2023 
'''
#%%
import numpy as np
import pandas as pd

TESS_id = np.loadtxt('observed_TIC-IDS.txt', dtype = str)

# --- Remove TIC from TESS ID ---
TESS_id = [id[3:] for id in TESS_id]
TESS_id = np.array(TESS_id, dtype = int)

# --- Loading TESS database ---
TESS_df = pd.read_csv('TESS_ext_target_data.csv')
TESS_df = TESS_df[TESS_df['TIC_ID'].isin(TESS_id)]

# --- Saving to .csv ---
TESS_df = TESS_df.to_csv('TESS_ext_target_data_observed.csv', index = False)