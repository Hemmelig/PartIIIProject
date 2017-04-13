#!/usr/bin/env python3

import obspy, os.path, glob, sys, scipy, scipy.signal
from obspy import read
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

event = sys.argv[1]

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')
print(seislist)

norm = None
azis = []

for s in range(len(seislist)):
    seis = read(seislist[s], format='PICKLE')
    sta  = seis[0].stats.station
    nw   = seis[0].stats.network
    slat = seis[0].stats['stla']
    slon = seis[0].stats['stlo']
    azis.append(seis[0].stats['az'])

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
    tshift = seis[0].stats['starttime'] -seis[0].stats['eventtime']
    
    Phase = ['Sdiff', 'S']
    for x in range(0, 2):
        if seis[0].stats.traveltimes[Phase[x]] != None:
            phase = Phase[x]
            if seis[0].stats['dist'] < 85:
                plt.subplot(1,4,1)
            elif seis[0].stats['dist'] < 90:
                plt.subplot(1,4,2)
            elif seis[0].stats['dist'] < 95:
                plt.subplot(1,4,3)
            elif seis[0].stats['dist'] < 100:
                plt.subplot(1,4,4)
                
    crs = open(station_name, "r")
    for columns in (row.strip().split() for row in crs):
        times.append(float(columns[0]) - 450 - seis[0].stats.traveltimes[phase])
        seistoplot.append(float(columns[1]))

    if norm == None:
        norm = 1. * np.max(seistoplot)
    seistoplot = [x / norm for x in seistoplot]
    plt.plot(times, seistoplot + np.round(seis[0].stats['az']), 'k', linewidth=1)

    plt.xlim([-20, 60])

    iInd = np.searchsorted(times, -20)
    fInd = np.searchsorted(times, 60)
    print(iInd, fInd, times[iInd], times[fInd])

    windowSeis = []
    windowTime = []
    for i in range(iInd, fInd):
        windowSeis.append(seistoplot[i])
        windowTime.append(times[i])

    analytical_signal = hilbert(windowSeis)
    amplitude_envelope = np.abs(analytical_signal)
    plt.plot(windowTime, amplitude_envelope + np.round(seis[0].stats['az']))

    # Finding peak maxima
    #peakind = scipy.signal.find_peaks_cwt(windowSeis, np.arange(1,100))
    #print(peakind, windowSeis[peakind[0]], windowTime[peakind[0]])

    #for x in range(len(peakind)):
    #    plt.plot(windowTime[peakind[x]], windowSeis[peakind[x]] + np.round(seis[0].stats['az']), 'x')
    
# Put labels on graphs
initDist = 85
for i in range(4):
    plt.ylim(round(min(azis) - 1), round(max(azis) + 1))
    plt.subplot(1,4,i + 1)
    plt.title('S dist < %d' % (initDist + (5 * i)), fontsize=10)
    if (i == 1):
        plt.ylabel('Azimuth (dg)')
    plt.xlabel('Time around predicted arrival (s)', fontsize=6)
    plt.tick_params(axis='both', which='major', labelsize=6)

# Save file and show plot
plt.savefig('Plots/' + event + '/' + event + '_synth_with_azimuth.pdf')

plt.show()
