'''
Code Purpose: Read in pickle files
Author: Owen A. Johnson
Last Modified: 2023-03-13 
'''
#%% 
import pickle
import numpy as np
import glob as glob
import pandas as pd

file_list = glob.glob('*.pkl')
freqs = np.arange(110, 200, 10)

# --- convert pickle object to pandas ---
Tsys_array = []; freq_array = []; ids = []

for file in file_list: 
    obj = pd.read_pickle(r'{}'.format(file))
    obj_df = pd.DataFrame(obj)
    keys = obj_df.keys()

    for key in keys:
        T_sys = obj_df[key][0][1]
    
        for f in freqs:
            Tsys_array.append(T_sys[f])
            freq_array.append(f)
            ids.append(key)

# --- Change column names ---
df = pd.DataFrame({'TICid':ids,'Tsys': Tsys_array, 'freq': freq_array})
df.to_csv('Tsys.csv', index=False)
