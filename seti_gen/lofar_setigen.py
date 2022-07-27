"""
Code Purpose: To inject signals into LOFAR data for the purposes of testing coincidence rejection. 
Author: Bryan Brzycki, edited by Owen Johnson.
Last Major Update: 27/07/2022
"""

from astropy import units as u
import setigen as stg
import blimpy as bl
import matplotlib.pyplot as plt

data_path = '/datax2/owen/bary.corrected.fils/TIC81831095.bary.0000.fil'
waterfall = bl.Waterfall(data_path, f_start = 145, f_stop = 155)
frame = stg.Frame(waterfall=waterfall)
frame.add_constant_signal(f_start=frame.get_frequency(150),
                          drift_rate=2*u.Hz/u.s,
                          level=frame.get_intensity(snr=30),
                          width=40*u.Hz,
                          f_profile_type='sinc2')

fig = plt.figure(figsize=(10, 6))
frame.plot()
plt.show()