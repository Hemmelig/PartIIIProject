#!/usr/bin/env python3

import obspy, obspy.signal, os.path, glob, scipy, sys
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

def findAlign(timeArray, seisArray):
    iInd = np.searchsorted(timeArray, -20)
    fInd = np.searchsorted(timeArray, 60)
    windowSeis = []
    for i in range(iInd, fInd):
        windowSeis.append(seisArray[i])
        
    maxInd = np.argmax(windowSeis) + iInd
    return maxInd

event = "20170224"
titles = ['a', 'b', 'c']

fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')

norm = None

fig = plt.figure(figsize=(6,6))

for s in range(1, 4, 1):
    seis = read(seislist[s], format='PICKLE')
    ax = fig.add_subplot(3,1,s)
    ax.set_title(titles[s - 1], fontsize=15, x=0.02, weight='bold')
    if (s == 2):
        sta = seis[0].stats.station
        nw = seis[0].stats.network
        slat = seis[0].stats['stla']
        slon = seis[0].stats['stlo']
        tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        while len(nw) < 2:
            nw = '_' + nw
        while len(sta) < 4:
            sta = '_' + sta
        if len(sta) > 4:
            sta = sta[:4]

        station_name = dir + 'Synthetics/UT_%s_%s' % (nw, sta)

        times = []
        seistoplot = []

        crs = open(station_name, 'r')
        for columns in (row.strip().split() for row in crs):
            times.append(float(columns[0]))
            seistoplot.append(float(columns[1]))

        norm = np.max(np.absolute(seistoplot))
        seistoplot = [x / norm for x in seistoplot]
        times = [x + tshift - seis[0].stats.traveltimes["S"] - 700 for x in times]

        alignInd = findAlign(times, seistoplot)

        talign = times[alignInd]

        times = [x - talign for x in times]

        ax.plot(times, seistoplot, 'k', linewidth=2)
       
    else:
        
        seistoplot = seis.select(channel='BHT')[0]
        seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)
        seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

        norm = 1. * np.max(abs(seistoplot.data))
        seistoplot.data = seistoplot.data / norm

        tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        seistoplot.times = [x + tshift - seis[0].stats.traveltimes["S"] for x in seistoplot.times()]

        alignInd = findAlign(seistoplot.times, seistoplot.data)

        talign = seistoplot.times[alignInd]
        seistoplot.times = [x - talign for x in seistoplot.times]

        ax.plot(seistoplot.times, seistoplot.data, 'k', linewidth=2)

    ax.set_ylim([-1.1, 1.1])
    ax.set_xlim([-15, 45])
    ax.set_ylabel('Normalised amplitude', fontsize=10)
    if (s == 3):
        ax.set_xlabel('Time around predicted arrival (s)', fontsize=10)
    else:
        ax.set_xticklabels([])

fig.savefig('Plots/overlap.pdf', bbox_inches='tight')
fig.savefig('Plots/overlap.eps', format='eps', dpi=1000, bbox_inches='tight')
        
plt.show()
