import numpy as np 
import matplotlib.pyplot as plt 
import csv
import pandas as pd 
import sys
import os
from matplotlib import colors

#vdfdThis code plots DM vs Time with SN being represented by colour

fig1=plt.figure(1)
ax1=fig1.add_subplot(111)

name=["P", "SN", "DM", "DMID", "NIDs", "fold_num", "P/Ptop", "Ptop/P", "file", "num"]
data=pd.read_csv("20180921.lis", sep='\s+', skiprows=0, names=name)
#index=data[data["S/N"] <= 9].index
#data=data.drop(index, inplace=False)
period=data["P"].values
DM=data["DM"].values
SN=data["SN"].values
print(DM)
pul=ax1.scatter(np.log10(period), DM, c=SN, cmap="viridis", vmin=np.min(SN), vmax=np.max(SN)) 
ax1.set_xlabel(" Log_10(Period) (ms)")
ax1.set_ylabel("DM (pc/cc)")
ax1.set_title("20180921 Periodicity")
ax1.set_xlim(2.3, 2.5)

cbar=fig1.colorbar(pul)
cbar.set_label("S/N")
plt.savefig("20180921_periodicity.png")


