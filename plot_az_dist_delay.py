#!/usr/bin/env python3

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import obspy
from obspy import read

# Event to plot
name = sys.argv[1]

cm = plt.cm.get_cmap('jet')

dir = 'Data/' + name + '/'
seislist = glob.glob(dir + '/*PICKLE')

delayTimes = []
amplRatios = []
azis = []
dists = []
dtdt = []

# Read in data from textfile
with open(dir + 'peakData.txt') as f:
    content = f.readlines()
    for x in content:
        var = x.split(':')
        delayTimes.append(float(var[0]))
        amplRatios.append(float(var[1]))
        azis.append(float(var[2]))
        dists.append(float(var[3]))

for s in range(len(seislist)):
    print(s + 1, len(seislist))
    seis = read(seislist[s], format='PICKLE')

    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    if (delayTimes[s] != 0):
        dtdt.append(delayTimes[s] - dtPred)
    else:
        dtdt.append(0)

# Make plot
#plot = plt.scatter(azis, dists, s=35, c=dtdt, vmin=10, vmax=16, alpha=1, cmap=cm)

#plt.colorbar(plot)

ds = np.arange(80, 101, 1)
azs = np.arange(5, 26, 1)
binAv = np.zeros((len(ds), len(azs)))
print(ds, azs)

for i in range(len(ds)):
    d = ds[i]
    for j in range(len(azs)):
        a = azs[j]
        for k in range(len(delayTimes)):
            avtmp = []
            if (azis[k] >= a and azis[k] <= a + 1 and dists[k] >= d and dists[k] <= d + 1):
                avtmp.append(dtdt[k])
            if (len(avtmp) > 0):
                binAv[i][j] = np.mean(avtmp)

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.set_title('Delay time as a function of distance and azimuth')
cax = ax.matshow(binAv, vmin=6, vmax=16)
fig.colorbar(cax)
#circle = plt.Circle()
#ax.add_artist(circle)

# Convert array values into strings
ds = [str(x) for x in ds]
azs = [str(x) for x in azs]

ax.set_xticklabels([''] + ds)
ax.set_yticklabels([''] + azs)

#plt.title('Delay time as a function of distance and azimuth')

plt.show()
