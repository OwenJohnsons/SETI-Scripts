
import numpy as np 
import matplotlib.pyplot as plt 
import csv
import pandas as pd 
import sys
import os
from matplotlib import colors

#This code plots DM vs Time with SN being represented by colour, using multiple pls files from multiple DM trials

fig1=plt.figure(1)
ax1=fig1.add_subplot(111)
DMs=np.linspace(0, 30, 1)


for file in os.listdir():
	filename = os.fsdecode(file)
	if filename.endswith(".pls"):
		names=["DM", "Duration", "time", "S/N", "length"]
		data=pd.read_csv(file, sep='\s+', skiprows=1, names=names) 
		DM=data["DM"].values
		time=(data["time"].values)*(163.84*(10**(-6)))
		SN=data["S/N"].values
		width=2**(data["Duration"].values)
		pul=ax1.scatter(time, DM,
		s=width, c=(SN), cmap="plasma",
		vmin=np.min(SN), vmax=np.max(SN))
		continue
	else:
		continue   
print(time)
ax1.set_xlabel("time (s)")
ax1.set_ylabel("DM (pc/cc)")
ax1.set_title(sys.argv[1])
cbar=fig1.colorbar(pul)
cbar.set_label("S/N")
plt.savefig("pic.png")
#plt.show()



