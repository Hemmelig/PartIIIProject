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

switch_yaxis = False

# Frequencies for filter
fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')

### Recorrect to original norm system once the data_selector is working
norm = None
dists = []
xlim = 0.

# Loop through seismograms
for s in range(len(seislist)):
    print(s + 1, len(seislist))
    seis = read(seislist[s], format='PICKLE')

    # List all distances
    dists.append(seis[0].stats['dist'])
    print(seis[0].stats['az'], seis[0].stats['dist'])
   
    seistoplot = seis.select(channel='BHT')[0]

    if (s == 0):
       xlim = seis[0].stats.traveltimes['S']

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
    seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
    seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)
    seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

    if syn:
        seissyn.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seissyn.filter('lowpass', freq=fmax, corners=2, zerophase=True)
        seissyn.data = np.gradient(seissyn.data, seissyn.stats.delta)

    # Time shift to shift data to reference time
    tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

    if norm == None:
        norm = 1. * np.max(abs(seistoplot.data))
    plt.plot(seistoplot.times() + tshift, seistoplot.data / norm + (seis[0].stats['dist']), 'k')

    if syn:
        plt.plot(seissyn.times(), seissyn.data / norm + (seis[0].stats['dist']), 'b')
      
    # Plot travel time predictions
    #for k in seis[0].stats.traveltimes.keys():
    #    if seis[0].stats.traveltimes[k] != None:
    #        plt.plot(seis[0].stats.traveltimes[k], np.round(seis[0].stats['dist']), 'g', marker='o', markersize=5)
    #        plt.text(seis[0].stats.traveltimes[k], np.round(seis[0].stats['dist'])-0.5, k, fontsize=6)  
   
# Put labels on graphs
plt.subplot(1,1,1)
plt.title(' ')
plt.ylabel('Distance (dg)')
plt.xlabel('Time around predicted arrival (s)')
plt.xlim([xlim - 100, xlim + 300])
plt.ylim([79, 98])
if switch_yaxis:
    plt.gca().invert_yaxis()

# Save file and show plot
if syn ==  True:
    plt.savefig('Plots/' + event + '/' + 'distance.pdf')
if syn == False:
    plt.savefig('Plots/' + event + '/' + 'distance_noSyn.pdf')

plt.show()
