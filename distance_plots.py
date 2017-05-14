#!/usr/bin/env python3

import obspy, obspy.signal, os.path, glob, scipy, sys
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

events = ['20170224', '20160528', '20160924']
titles = ['a', 'b', 'c']

switch_yaxis = False

def findAlign(timeArray, seisArray):
    iInd = np.searchsorted(timeArray, -20)
    fInd = np.searchsorted(timeArray, 60)
    windowSeis = []
    for i in range(iInd, fInd):
        windowSeis.append(seisArray[i])

    maxInd = np.argmax(windowSeis) + iInd
    return maxInd

fig = plt.figure(figsize=(12,12))

for i in range(len(events)):
    # Frequencies for filter
    if (i == 2):
        fmin = 0.033
        fmax = 0.05
    else:
        fmin = 0.033
        fmax = 0.1

    dir = 'Data/' + events[i] + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    norm = None

    dists = []
    
    ax = fig.add_subplot(1,4,i+1)

    for s in range(len(seislist)):
        print(s + 1, len(seislist))
        seis = read(seislist[s], format='PICKLE')

        # List all distances
        dists.append(seis[0].stats['dist'])
   
        seistoplot = seis.select(channel='BHT')[0]

        # Plot seismograms
        phase = 'S'

        # Filter data
        seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)
        seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

        # Time shift to shift data to reference time
        tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        seistoplot.times = [x + tshift - seis[0].stats.traveltimes['S'] for x in seistoplot.times()]

        alignInd = findAlign(seistoplot.times, seistoplot.data)

        talign = seistoplot.times[alignInd]

        seistoplot.times = [x - talign for x in seistoplot.times]
    
        if norm == None:
            if (i < 2):
                norm = 0.8 * np.max(abs(seistoplot.data))
            else:
                norm = 1. * np.max(abs(seistoplot.data))
        ax.plot(seistoplot.times, seistoplot.data / norm + np.round((seis[0].stats['dist'])), 'k')

        # Put labels on graphs
        if (i == 0):
            ax.set_ylabel('Distance (dg)', fontsize=16)
        else:
            ax.set_yticklabels([])
        ax.set_xlim([-15, 50])
        ax.set_ylim([78, 100])
        ax.set_title(titles[i], x=0.1, weight='bold', fontsize=14)


dir = 'Data/' + events[0] + '/'
seislist = glob.glob(dir + '/*PICKLE')

norm = None
dists = []
    
ax = fig.add_subplot(1,4,4)

for s in range(len(seislist)):
    seis = read(seislist[s], format='PICKLE')
    sta  = seis[0].stats.station
    nw   = seis[0].stats.network
    slat = seis[0].stats['stla']
    slon = seis[0].stats['stlo']
    dists.append(seis[0].stats['dist'])

    while len(nw) < 2:
        nw = '_' + nw
    while len(sta) < 4:
        sta = '_' + sta
    if len(sta) > 4:
        sta = sta[:4]

    station_name = dir + 'Synthetics/UT_%s_%s' % (nw, sta)
    print(station_name)

    times = []
    seistoplot = []

    # Time shift to shift data to reference time
    tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']
    
    phase = 'S'
    
    crs = open(station_name, "r")
    for columns in (row.strip().split() for row in crs):
        times.append(float(columns[0]) - 450 - seis[0].stats.traveltimes[phase])
        seistoplot.append(float(columns[1]))

    if norm == None:
        norm = 0.6 * np.max(seistoplot)
    seistoplot = [x / norm for x in seistoplot]
    ax.plot(times, seistoplot + np.round(seis[0].stats['dist']), 'k', linewidth=1)

    # Put labels on graphs
    ax.set_yticklabels([])
    ax.set_xlim([-15, 50])
    ax.set_ylim([78, 100])
    ax.set_title('d', x=0.1, weight='bold', fontsize=14)

fig.text(0.5, 0.04, 'Time around predicted arrival (s)', ha='center', fontsize=16)    
    
plt.savefig('Plots/distances.pdf', bbox_inches='tight')
plt.savefig('Plots/distances.eps', format='eps', dpi=1000, bbox_inches='tight')

plt.show()
