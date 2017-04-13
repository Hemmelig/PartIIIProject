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

# Frequencies for filter
fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')

### Recorrect to original norm system once the data_selector is working
normReal, normSyn = None, None
dists = []

# Loop through seismograms
for s in range(0, len(seislist), 4):
    print(s, len(seislist))
    seis = read(seislist[s], format='PICKLE')

    # List all distances
    dists.append(seis[0].stats['dist'])
    print(seis[0].stats['az'], seis[0].stats['dist'])
   
    seistoplot = seis.select(channel='BHT')[0]

    # Plot synthetics
    if syn:
        seissyn = seis.select(channel='BXT')[0]

    # Plot seismograms
    print(seis[0].stats.traveltimes['Sdiff'])
    Phase = ['Sdiff', 'S']
    for x in range(0,2):      
        if seis[0].stats.traveltimes[Phase[x]] != None:
            phase = Phase[x]
            plt.subplot(1,1,1)

    # Filter data
    if real:
        seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)
    if syn:
        seissyn.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seissyn.filter('lowpass', freq=fmax, corners=2, zerophase=True)

    # Time shift to shift data to reference time
    tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']
    print('max', np.max(seistoplot.times()))

    if real:
        if normReal == None:
            normReal = -1. * np.max(seistoplot.data)
        plt.plot(seistoplot.times() + tshift, seistoplot.data / normReal + (seis[0].stats['dist']), 'k')

    if syn:
        if normSyn == None:
            #normSyn = 1. * np.max(seissyn.data)
            normSyn = 1. * np.max(seistoplot.data)
        print(np.max(seistoplot.data), np.max(seissyn.data))
        plt.plot(seissyn.times(), seissyn.data / normSyn + (seis[0].stats['dist']), 'b')
      
    # Plot travel time predictions
    for k in seis[0].stats.traveltimes.keys():
        if seis[0].stats.traveltimes[k] != None:
            plt.plot(seis[0].stats.traveltimes[k], np.round(seis[0].stats['dist']), 'g', marker='o', markersize=5)
            plt.text(seis[0].stats.traveltimes[k], np.round(seis[0].stats['dist'])-0.5, k, fontsize=6)
   
# Put labels on graphs
plt.subplot(1,1,1)
plt.title(' ')
plt.ylabel('Distance (dg)')
plt.xlabel('Time around predicted arrival (s)')
plt.xlim([1250, 1400])
plt.ylim([70, 100])
if switch_yaxis:
    plt.gca().invert_yaxis()

# Save file and show plot
if syn ==  True:
    plt.savefig('Plots/' + event + '/' + event + '_data_with_distance.pdf')
if syn == False:
    plt.savefig('Plots/' + event + '/' + event + '_data_with_distance_noSyn.pdf')

plt.show()
