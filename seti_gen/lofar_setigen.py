"""
Code Purpose: To inject signals into LOFAR data for the purposes of testing coincidence rejection. 
Author: Bryan Brzycki, edited by Owen Johnson.
Last Major Update: 27/07/2022
"""

from astropy import units as u
import setigen as stg
import blimpy as bl
import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family = 'serif', serif = 'cmr10') 
plt.rcParams.update({'font.size': 22})
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'stix'

fstr = 153.00025; fstop = 153.00075 # - aribitrary frequency width of 250 Hz for testing purposes
mid_f = round(np.abs(fstr + fstop)/2.,  4)

# --- Loading data files and loading to a frame for injection --- 
print('Loading data....')
data_path = '/datax2/owen/bary.corrected.fils/TIC81831095.bary.0000.fil'
waterfall = bl.Waterfall(data_path, f_start = fstr, f_stop = fstop) # - Taking a frquency subset 

data_shape = waterfall.data.shape
frame = stg.Frame(waterfall=waterfall)

print('Starting signal injection...')


# --- Actual Signal Injection ---

drift_rates = np.random.uniform(-6,6,[10])
sn_values = np.random.uniform(50, 100,[10])

# dr = 2.711
# sr = 88

for i in range(0, len(drift_rates)):
    frame.add_constant_signal(f_start=frame.get_frequency(data_shape[2]/3),
                            drift_rate=drift_rates[i]*u.Hz/u.s,
                            level=frame.get_intensity(snr=sn_values[i]),
                            width=3*u.Hz,
                            f_profile_type='sinc2')


plt.style.use('dark_background')

fig = plt.figure(figsize=(10, 8), dpi = 300)
frame.bl_plot()

plt.title(('TIC81831095 \n  Injected Signal (SE)'), fontsize = 24)
plt.xlabel('Relative Frequency [Hz] from %s [MHz]' % mid_f); plt.ylabel('Time [s]')

# --- x-axis management --- 
factor = 1e6 # - Scale for plotting purposes
xloc = np.linspace(fstr, fstop, 5)
xticks = [round(loc_freq) for loc_freq in (xloc - mid_f)*factor]
if np.max(xticks) > 1000:
	xticks = [xt/1000 for xt in xticks]
	units = 'Hz'
plt.xticks(xloc, xticks)
plt.savefig('/datax2/owen/seti_gen.png', bbox_inches='tight', transparent = True)
plt.show()

frame.save_fil('generated_.h5')

print('Injecting Done!')