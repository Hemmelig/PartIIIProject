#!/usr/bin/env python3

import sys, obspy, obspy.signal, os.path, glob, scipy, subprocess
from obspy import read
from obspy.core import Stream
import numpy as np

heights = np.arange(35, 3, -2)

def timeToHeight(dt, dist):
    with open(dir + '_' + str(int(dist)) + '_2d_dt.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        dt_arr = [[float(e) for e in r] for r in reader]

    tmp_arr = []
    for i in range(len(dt_arr[0][0])):
        tmp_arr.append(dt_arr[i][5])
    idx = np.abs(tmp_arr - dt).argmin()

    del tmp_arr
    
    return heights[idx]

dir = 'Data/' + name + '/'
seislist = glob.glob(dir + '/*PICKLE')

# Read in data from textfile
f = open(dir + 'peakData.txt', 'r')

delayTimes = []
dists = []
convHeights = []

pCorrection = 14.6 / 2.

lines = f.readlines()
for x in lines:
    delayTimes.append(float(x.split(':')[0]))
    dists.append(float(x.split(':')[3]))
f.close()

# Loop through stations
for s in range(len(seislist)):
    print(s + 1, len(seislist[s]), seislist[s])
    seis = read(seislist[s], format='PICKLE')

    dist = np.round(dists[s])
    
    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    
    # Get event-station pair delay time and amplitude ratio
    delayTime = delayTimes[s]
    if (delayTime == 0):
        convHeights.append(0)
    else:
        delayTime = delayTime - dtPred - pCorrection
        height = timeToHeight(delayTime, dist)
        convHeights.append(height)

heightData = open(dir + 'convHeights.txt', 'w')
for i in range(len(convHeights)):
    heightData.write("%s \n" % (convHeights[i]))
heightData.close()

