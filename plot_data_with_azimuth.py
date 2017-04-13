#!/usr/bin/env python3

import obspy, obspy.signal, os.path, glob, scipy, sys
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

event = sys.argv[1]

synthetics = ''
while (synthetics != 'y' and synthetics != 'n'):
    synthetics = input('Do you want to plot synthetics (y/n)? ')
    if (synthetics == 'y'):
        syn = True
    if (synthetics == 'n'):
        syn = False

real = True # Plot real data
switch_yaxis = False

# Frequencies for filter in Hz
fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')

norm = None
azis = []

# Loop through seismograms
for s in range(0, len(seislist), 4):
    print(s, len(seislist))
    seis = read(seislist[s], format='PICKLE')

    # List all azimuths
    azis.append(seis[0].stats['az'])
    print(seis[0].stats['az'], seis[0].stats['dist'])

    # Plot synthetics
    if syn:
        seissyn = seis.select(channel='BXT')[0]

    # Split seismograms by distance range
    print(seis[0].stats.traveltimes['Sdiff'])
    Phase = ['Sdiff', 'S']
    for x in range(0,2):
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

    # Filter data
    if real:
        seis.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seis.filter('lowpass', freq=fmax, corners=2, zerophase=True)
    if syn:
        seissyn.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seissyn.filter('lowpass', freq=fmax, corners=2, zerophase=True)

    seistoplot = seis.select(channel='BHT')[0]
    scalerD = np.max(seistoplot.data)

    # Differentiate data streams to produce velocity seismograms
    seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)
    if syn:
        seissyn.data = np.gradient(seissyn.data, seissyn.stats.delta)

    # Time shift to shift data to reference time
    tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']
    print('max', np.max(seistoplot.times()))

    if real:
        if norm == None:
            norm = 0.25 * scalerD
        plt.plot(seistoplot.times() + tshift - seis[0].stats.traveltimes[phase], seistoplot.data / norm + np.round(seis[0].stats['az']), 'k', linewidth=1)
    if syn:
        if norm == None:
            norm = 0.25 * scalerD
        plt.plot(seissyn.times() - seis[0].stats.traveltimes[phase], seissyn.data / norm + np.round(seis[0].stats['az']), 'b', linewidth=1)

    plt.xlim([-20, 60])

    # Plot travel time predictions
    for k in seis[0].stats.traveltimes.keys():
        if seis[0].stats.traveltimes[k] != None:
            plt.plot(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az']), 'g', marker='o', markersize=4)
            # Uncomment line below if wish to display phase names ***TO IMPLEMENT*** Question at start to question if wanted
            #plt.text(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az'])-0.5, k, fontsize=6)           
            
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
    if switch_yaxis:
        plt.gca().invert_yaxis()
   
# Save file and show plot
if syn ==  True:
    plt.savefig('Plots/' + event + '/' + event + '_data_with_azimuth_v.pdf')
if syn == False:
    plt.savefig('Plots/' + event + '/' + event + '_data_with_azimuth_noSyn_v.pdf')

plt.show()     

