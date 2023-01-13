'''
Code Purpose: Intake a text file and query the gaia database for targets with matching co-ordinates and returing a .csv for the purposes of further analysis. 
Author: Owen A. Johnson
Last Major Update: 11/01/2023
'''

import argparse
import sqlite3
import time
import numpy as np 
import pandas as pd
from tqdm import tqdm

start = time.time()

# --- Arguments ---
parser = argparse.ArgumentParser(description='Intake a text file and query the gaia database for targets with matching Gaia database ids and returing a .csv for the purposes of further analysis.')
parser.add_argument('-i', '--inputfile', type=str, help='Text file gaia source ids to be indexed.')

# --- Loading targets from input file --- 
input_ids = np.loadtxt(str(parser.parse_args().inputfile), unpack = True) # Unpacking targets found from Gaia in observed LOFAR pointings in RA and DEC
print('Number of targets to be indexed: ', len(input_ids))

# --- Connecting to database ---
connection = sqlite3.connect('gaia_master_database.db')
connection.row_factory = sqlite3.Row
cursor = connection.execute('select * from data') # - fetching 'data' table from database
# cursor = connection.execute("SELECT * FROM data LIMIT 100000" ) # limits the query to 'n' entries; for testing purposes. 

# --- fetching column names ---
row = cursor.fetchone()
rows = cursor.fetchall()
names = row.keys()
print(names)


# --- Creating a dataframe of the matches ---
dataframe = pd.read_sql_query(f'SELECT * FROM data WHERE "Source_id" IN {tuple(input_ids)}', connection) # Create a DataFrame from the table
dataframe = dataframe.drop(columns = ['index']) # Drop the index column as it is not needed
print(dataframe.head())
dataframe.to_csv('gaia_inbeam_table.csv',index = False)

# Close the connection
cursor.close()
end = time.time()

print('Time taken: %s seconds. ' % '{:.2f}'.format(end - start))