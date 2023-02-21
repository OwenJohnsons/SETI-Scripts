'''
Title: Query Vizier with unique Gaia IDs and return a dataframe. 
Code Purpose: 
Last Major Update: 21/02/2022 
'''

import pandas as pd 
import numpy as np 
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
Vizier.ROW_LIMIT = -1

beam_ids = np.loadtxt('txt/gaia_ids.txt')

# --- Using astroquery.gaia to query with the array of IDs ---
query = """SELECT * FROM gaiadr2.gaia_source WHERE source_id IN %s""" % (tuple(beam_ids),)
job = Gaia.launch_job(query)
print(job)