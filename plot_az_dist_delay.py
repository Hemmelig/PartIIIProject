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
bounceLats = []
bounceLons = []
pCorrection = 29.13 / 4.

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

    turnloc = seis[0].stats.turnpoints["ScS"]
    bounceLats.append(turnloc[1])
    bounceLons.append(turnloc[2] + 360.)

    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    if (delayTimes[s] != 0):
        dtdt.append(delayTimes[s] - dtPred - pCorrection)
    else:
        dtdt.append(0)

lons = np.arange(184, 204, 2)
lats = np.arange(8, 28, 2)
binAv = np.zeros((len(lons), len(lats)))
binAv.fill(-1)

for i in range(len(lons)):
    lon = lons[i]
    for j in range(len(lats)):
        lat = lats[j]
        for k in range(len(delayTimes)):
            avtmp = []
            if (bounceLats[k] >= lat and bounceLats[k] <= lat + 2 and bounceLons[k] >= lon and bounceLons[k] <= lon + 2):
                avtmp.append(dtdt[k])
            if (len(avtmp) > 0):
                binAv[j][i] = np.mean(avtmp)

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

#print(binAv)

LATS, LONS = np.meshgrid(lats, lons)

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.set_title('Delay time map ' + name)
binAv[np.isnan(binAv)] = -1
binAv = np.ma.array(binAv, mask=binAv == -1)
plt.gca().patch.set_color('.25')
cax = plt.contourf(binAv)
#cax = plt.contourf(LONS, LATS, binAv)
circle = plt.Circle((192.5, 17.5), 7.5, color="white", linewidth=2, fill=False)
plt.gca().add_artist(circle)
#cax = ax.matshow(binAv, vmin=2, vmax=8)
fig.colorbar(cax)

# Convert array values into strings
los = [str(x) for x in lons]
las = [str(x) for x in lats]

ax.set_xticklabels(los)
ax.set_yticklabels(las)

plt.show()
