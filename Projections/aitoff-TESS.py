''' 
Code Purpose: Aitoff projection of TESS pointings with LOFAR beam size
Author: Owen A. Johnson 
Last Major Update: 2023-03-27
'''
#%% 
import numpy as np
import matplotlib.pyplot as plt
import scienceplots; plt.style.use('science')
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

color_pallette = ['#016f95', '#0098c0', '#28c062', '#f05039', '#bd2d1d']

ext_df = pd.read_csv('csv/TESS_ext_target_data.csv'); obs_df =  pd.read_csv('csv/observed.targets.csv')
dropped_df = ext_df[~ext_df.TIC_ID.isin(obs_df['TIC_ID'])] # Dropping entries that have been observed to avoid duplicates

# - Loading Co-ordinates into Arrays - 
obs_ra = np.array(obs_df['ra']); obs_dec = np.array(obs_df['dec'])
all_ra = np.array(dropped_df['ra']); all_dec = np.array(dropped_df['dec'])

obs_trgts = SkyCoord(obs_ra, obs_dec, frame='fk5', unit=(u.hourangle, u.deg))
all_trgts = SkyCoord(all_ra, all_dec, frame = 'fk5', unit=(u.hourangle, u.deg))

gal_obs = obs_trgts.galactic
gal_all = all_trgts.galactic

# Convert the circle points to Aitoff coordinates
npoints = 100
theta = np.linspace(0, 2*np.pi, npoints)
x = np.cos(theta); y = np.sin(theta)
r = 1.295 # Radius of LOFAR beam in degrees

fig = plt.figure(figsize=(10, 5), dpi=200)
ax = fig.add_subplot(111, projection='aitoff')

for i in range(0, len(gal_obs)):
    l = gal_obs[i].l.deg; b = gal_obs[i].b.deg
    gal_points = SkyCoord(l + r*x/np.cos(np.radians(b)), b + r*y, unit='deg', frame='galactic')
    ax.plot(gal_points.l.wrap_at('180d').radian, gal_points.b.radian, color=color_pallette[2], linewidth=0.2, zorder = 5)
    
ax.plot(gal_points.l.wrap_at('180d').radian, gal_points.b.radian, color=color_pallette[2], linewidth=0.2, zorder = 5, label = 'LOFAR Beam')

ax.scatter(gal_obs.l.wrap_at('180d').radian, gal_obs.b.radian, color=color_pallette[3], s=0.2, label = 'Observed Pointings', zorder = 5)
ax.scatter(gal_all.l.wrap_at('180d').radian, gal_all.b.radian, facecolor = 'None', edgecolors='grey', s = 5, label = 'Targets of Interest', zorder = 2, alpha = 0.2)

plt.axhspan(np.deg2rad(-5), np.deg2rad(5), alpha=0.2, color='grey', label = 'Galactic Center')

plt.ylabel('Galactic Latitude', fontsize = '18'); plt.xlabel('Galactic Longitude', fontsize = '18')
ax.yaxis.tick_right()
plt.legend(frameon = True)
plt.grid()
plt.savefig('Aitoff_Projection.png', transparent=True, bbox_inches='tight', dpi = 500)
plt.show()