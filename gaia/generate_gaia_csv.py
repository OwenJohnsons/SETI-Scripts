'''
Code Purpose: Intake a text file and query the gaia database for targets with matching co-ordinates and returing a .csv for the purposes of further analysis. 
Author: Owen A. Johnson
Last Major Update: 11/01/2023
'''

import argparse
import sqlite3
import numpy as np 
import pandas as pd
from tqdm import tqdm

# --- Arguments ---
parser = argparse.ArgumentParser(description='Intake a text file and query the gaia database for targets with matching co-ordinates and returing a .csv for the purposes of further analysis.')
parser.add_argument('-i', '--inputfile', type=str, help='Text file containing RA and DEC co-ordinates of targets to be queried from the gaia database.')

# --- Loading targets from input file --- 
input_ra, input_dec = np.loadtxt(str(parser.parse_args().inputfile), unpack = True) # Unpacking targets found from Gaia in observed LOFAR pointings in RA and DEC

# --- Connecting to database ---
connection = sqlite3.connect('gaia_master_database.db')
connection.row_factory = sqlite3.Row
cursor = connection.execute('select * from data') # - fetching 'data' table from database
# cursor = connection.execute("SELECT * FROM data LIMIT 1000" ) # limits the query to 'n' entries; for testing purposes. 

# --- fetching column names ---
row = cursor.fetchone()
rows = cursor.fetchall()
names = row.keys()
print(len(names))

# --- Finding RA and DEC matches in the database --- 

# The index of the column to check
column_index = 2

# Iterate through the rows
indexes = []
for i, row in tqdm(enumerate(rows), leave = True):
    if row[4] in input_ra and row[6] in input_dec: # 4 and 6 are the column indexes for ra and dec but the names can be used respectivly if preferred. 
        indexes.append(i)


# --- Creating a dataframe of the matches ---

dataframe = pd.read_sql_query(f'SELECT * FROM data WHERE "index" IN {tuple(indexes)}', connection) # Create a DataFrame from the table
dataframe = dataframe.drop(columns = ['index']) # Drop the index column 
dataframe = dataframe.drop_duplicates() # Drop any duplicate entries
print(dataframe.head())
dataframe.to_csv('gaia_inbeam_table.csv',index = False)

# Close the connection
cursor.close()