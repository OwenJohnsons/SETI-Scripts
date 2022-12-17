'''
Code Purpose: To turn large .csv into sqlite database
Author: Owen A. Johnson 
Last Major Update: 16/12/2022
'''

import pandas as pd
import sqlite3
import time
import argparse
from sqlalchemy import create_engine
from tqdm import tqdm 

# - Command line arguments -
parser = argparse.ArgumentParser(description='Turn large .csv into sqlite database')
parser.add_argument('-i', '--input', type=str, help='Input .csv file', required=True)
parser.add_argument('-o', '--output', type=str, help='Output .db file', required=True)
args = parser.parse_args()

# - Loading .csv file -
print(pd.read_csv(args.input, nrows=5)) # Printing first 5 rows of .csv file
database = create_engine('sqlite:///' + args.output + '.db') # Creating database

# - looping over chunks of a manageable size and appending to database -
i = 0
j = 0 
for chunk in tqdm(pd.read_csv(args.input, chunksize=1e5, iterator = True)): # Reading in chunks of 100,000 rows at a time
    chunk = chunk.rename(columns={c: c.replace(' ', '') for c in chunk.columns}) # Removing spaces from column names
    chunk.index += j; i += 1
    chunk.to_sql('data', database, if_exists='append') # Appending to database
    j = chunk.index[-1] + 1
