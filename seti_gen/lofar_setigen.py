"""
Code Purpose: To inject signals into LOFAR data for the purposes of testing coincidence rejection. 
Author: Bryan Brzycki, edited by Owen Johnson.
Last Major Update: 27/07/2022
"""

from astropy import units as u
import setigen as stg
import blimpy as bl
import matplotlib.pyplot as plt

plt.rc('font', family = 'serif', serif = 'cmr10') 
plt.rcParams.update({'font.size': 22})

data_path = '/datax2/owen/bary.corrected.fils/TIC81831095.bary.0000.fil'
waterfall = bl.Waterfall(data_path, f_start = 153.00025, f_stop = 153.00050)

fig = plt.figure(figsize=(10, 8), dpi = 300)
waterfall.plot_spectrum(logged=True)
plt.savefig('raw_spectrum.png')
plt.show()

fig = plt.figure(figsize=(10, 8), dpi = 300)
waterfall.plot_waterfall()
plt.title('TIC81831095 (RAW)')

plt.savefig('raw_waterfall.png')
plt.show()

print('Starting signal injection...')

data_shape = waterfall.data.shape


frame = stg.Frame(waterfall=waterfall)
frame.add_constant_signal(f_start=frame.get_frequency(data_shape[2]/3),
                          drift_rate=0.2*u.Hz/u.s,
                          level=frame.get_intensity(snr=250),
                          width=5*u.Hz,
                          f_profile_type='sinc2')

fig = plt.figure(figsize=(10, 8), dpi = 300)
frame.plot()

# plt.title('TIC81831095, Injected (IE)', fontsize = 24)
plt.xlabel('Frequency [Hz]'); plt.ylabel('Time [s]')
plt.savefig('seti_gen.png', bbox_inches='tight')
plt.show()

print('Injecting Done!')