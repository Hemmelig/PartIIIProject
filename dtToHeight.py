#!/usr/bin/env python3

import sys, obspy, obspy.signal, os.path, glob, scipy, subprocess
from obspy import read
from obspy.core import Stream
import numpy as np
import csv

heights = np.arange(45, 3, -1)

def timeToHeight(dt, dist, hs):
    dir1 = 'Data/PeakData/20170224'
    if (dist > 97):
        dist = 97
        
    with open(dir1 + '_' + str(dist) + '_2d_dt.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        dt_arr = np.array([[float(e) for e in r] for r in reader])
            
    tmp_arr = np.array(dt_arr[:,5])
    print(tmp_arr)
    if (dt > max(tmp_arr)):
        dt = max(tmp_arr)
    if (dt < 3):
        return 0
        
    idx = np.argmin(np.abs(tmp_arr - dt))
    
    del tmp_arr
    
    return hs[idx]

names = ['20170224', '20160528']
pCorrections = [7.3, 7.6]

convHeights = []
tlats = []
tlons = []

for i in range(len(names)):
    name = names[i]
    pCorrection = pCorrections[i]

    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    # Read in data from textfile
    f = open(dir + 'peakData.txt', 'r')

    delayTimes = []
    dists = []
    pCorrection = 14.6 / 2.

    lines = f.readlines()
    for x in lines:
        delayTimes.append(float(x.split(':')[0]))
        dists.append(float(x.split(':')[3]))
    f.close()

    # Loop through stations
    for s in range(len(seislist)):
        print(s + 1, len(seislist), seislist[s])
        seis = read(seislist[s], format='PICKLE')

        turnloc = seis[0].stats.turnpoints['ScS']

        tlats.append(turnloc[1])
        tlons.append(turnloc[2])

        dist = np.round(dists[s])

        if (i == 0):
            dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
        else:
            dtPred = 0.
    
        # Get event-station pair delay time and amplitude ratio
        delayTime = delayTimes[s]
        if (delayTime == 0):
            convHeights.append(0)
        else:
            delayTime = delayTime - dtPred - pCorrection
            print(delayTime)
            height = timeToHeight(delayTime, int(dist), heights)
            convHeights.append(height)

dir = 'Data/'
        
heightData = open(dir + 'convHeights_20.txt', 'w')
for i in range(len(convHeights)):
    heightData.write("%s:%s:%s \n" % (convHeights[i], tlats[i], tlons[i]))
    print('Printing: ' + str(i))
heightData.close()
print('Done')
