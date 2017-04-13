#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy, obspy.signal, os.path, time, glob, shutil, scipy, sys, subprocess
from obspy import UTCDateTime, read
from obspy.core import Stream, Trace, event
from obspy.taup.taup import getTravelTimes
import matplotlib.pyplot as plt
import numpy as np
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call
import matplotlib.colors as colors
import matplotlib.cm as cm

#####################################################################################

def plotWithDistance(name):
    # Settings - may want to include this in the GUI once this script is done, but
    # as it is, each event will require significant editing anyway
    syn = True
    real = True
    colour = False
    switch_yaxis = False

    # Frequencies for filter in Hz
    fmin = 0.033
    fmax = 0.1

    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    savedir = 'Plots/' + name + '/'
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    
    norm = None
    dists = []

    # Loop through seismograms
    for s in range(0, len(seislist), 4):
        try:
            print(s)
            seis = read(seislist[s], format='PICKLE')
            dists.append(seis[0].stats['dist'])
            print(seis[0].stats['az'], seis[0].stats['dist'])
            seistoplot = seis.select(channel='BHT')[0]

            # Plot synthetics
            if syn:
                seissyn = seis.select(channel='BXT')[0]

            # Split seismograms by distance range (Needs to be adapted per event to produce a reasonable plot
            print(seis[0].stats.traveltimes['Sdiff'])
            Phase = ['Sdiff', 'S']

            for x in range(0,1):      
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
                if norm == None:
                    norm = 1. * np.max(seistoplot.data)
                plt.plot(seistoplot.times() + tshift, seistoplot.data / norm + np.round(seis[0].stats['dist']), 'k')

            if syn:
                if norm == None:
                    norm = 1. * np.max(seissyn.data)
                plt.plot(seissyn.times(), seissyn.data / norm + np.round(seis[0].stats['dist']), 'b')
            
            print(np.max(seistoplot.data), np.max(seissyn.data))
            #plt.xlim([-30, 40])
            # Plot travel time predictions
            for k in seis[0].stats.traveltimes.keys():
                if seis[0].stats.traveltimes[k] != None:
                    plt.plot(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az']), 'g', marker='o', markersize=4)

        except:
            pass
                
    # Put labels on graphs
    plt.subplot(1, 1, 1)
    plt.title(' ')
    plt.ylabel('Distance (dg)')
    plt.xlabel('Time around predicted arrival (s)')
    plt.xlim([0, 2700])
    plt.ylim([60, 110])
    if switch_yaxis:
        plt.gca().invert_yaxis()      

    # Save file and show plot
    #plt.savefig('Plots/' + name + '/' + name + '_data_with_distance.pdf')
    plt.show()
