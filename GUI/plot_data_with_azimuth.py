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

def plotWithAzimuth(name):
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
    azis = []

    # Plot measured differential travel times by colouring the lines
    if colour:
        dtsta = []
        dtnw = []
        dtdt = []
        rd = open('traveltimes_xcorr_12_30s_20130715.dat', 'r')
        for line in rd.readlines():
            val = line.split()
            dtsta.append(val[0])
            dtnw.append(val[1])
            dtdt.append(float(val[4]))
            cNorm = colors.Normalize(vmin=-10, vmax=10)
            scalarMap = cm.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet'))

    # Loop through seismograms
    for s in range(0, len(seislist), 4):
      #  try:
            print(s)
            seis = read(seislist[s], format='PICKLE')
            azis.append(seis[0].stats['az'])
            print(seis[0].stats['az'], seis[0].stats['dist'])
            seistoplot = seis.select(channel='BHT')[0]

            # Plot synthetics
            if syn:
                seissyn = seis.select(channel='BXT')[0]

            # Split seismograms by distance range (Needs to be adapted per event to produce a reasonable plot
            print(seis[0].stats.traveltimes['Sdiff'])
            Phase = ['Sdiff', 'S']
            
            for x in range(0, 2):
                if seis[0].stats.traveltimes[Phase[x]] != None:
                    phase = Phase[x]
                    if seis[0].stats['dist'] < 85:
                        plt.subplot(1, 4, 1)
                    elif seis[0].stats['dist'] < 90:
                        plt.subplot(1, 4, 2)
                    elif seis[0].stats['dist'] < 95:
                        plt.subplot(1, 4, 3)
                    elif seis[0].stats['dist'] < 100:
                        plt.subplot(1, 4, 4)

            # Filter data
   
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
                    norm = 2. * np.max(seistoplot.data)
                plt.plot(seistoplot.times() + tshift - seis[0].stats.traveltimes[phase], seistoplot.data / norm + np.round(seis[0].stats['az']), 'k')

            if syn:
                if norm == None:
                    norm = 2. * np.max(seissyn.data)
                plt.plot(seissyn.times() - seis[0].stats.traveltimes[phase], seissyn.data / norm + np.round(seis[0].stats['az']), 'b')
            print(np.max(seistoplot.data), np.max(seissyn.data))
            plt.xlim([-30, 40])
            # Plot travel time predictions
            for k in seis[0].stats.traveltimes.keys():
                if seis[0].stats.traveltimes[k] != None:
                    plt.plot(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az']), 'g', marker='o', markersize=4)

     #   except:
     #       pass
                    
    # Put labels on graphs
    plt.subplot(1, 4, 1)
    plt.title('S dist < 90')
    plt.ylabel('azimuth (dg)')
    if switch_yaxis:
        plt.gca().invert_yaxis()

    plt.ylim(round(min(azis) - 1), round(max(azis) + 1))
    plt.subplot(1, 4, 2)
    plt.title('Sdiff dist < 100')
    plt.xlabel('Time around predicted arrival (s)')
    if switch_yaxis:
        plt.gca().invert_yaxis()

    plt.ylim(round(min(azis) - 1), round(max(azis) + 1))
    plt.subplot(1, 4, 3)
    plt.title('Sdiff dist < 115')
    plt.xlabel('Time around predicted arrival (s)')
    if switch_yaxis:
        plt.gca().invert_yaxis()

    plt.ylim(round(min(azis) - 1), round(max(azis) + 1))
    plt.subplot(1, 4, 4)
    plt.title('Sdiff dist < 120')
    plt.xlabel('Time around predicted arrival (s)')
    if switch_yaxis:
        plt.gca().invert_yaxis()       

    # Save file and show plot
    plt.savefig('Plots/' + name + '/' + name + '_data_with_azimuth.pdf')
    plt.show()
