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

def plotWithAzimuth(name, syn):
    colour = False
    switch_yaxis = False

    # Frequencies for filter in Hz
    fmin = 0.033
    fmax = 0.1

    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')

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
    for s in range(len(seislist)):
        print(s, len(seislist))

        # Retrieve, filter and differentiate data
        seis = read(seislist[s], format='PICKLE')
        seis.filter('highpass', freq=fmin, corners=2, zerophase=True)
        seis.filter('lowpass', freq=fmax, corners=2, zerophase=True)
        seistoplot = seis.select(channel='BHT')[0]
        seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

        # list all azimuths
        azis.append(seis[0].stats['az'])

        # Split seismograms by distance range
        phase = 'S'
        if seis[0].stats['dist'] < 85:
            plt.subplot(1,4,1)
        elif seis[0].stats['dist'] < 90:
            plt.subplot(1,4,2)
        elif seis[0].stats['dist'] < 95:
            plt.subplot(1,4,3)
        elif seis[0].stats['dist'] < 100:
            plt.subplot(1,4,4)

        # Normalise data and shift data up by azimuth
        norm = 1. * np.max(abs(seistoplot.data))
        seistoplot.data = seistoplot.data / norm + np.round(seis[0].stats['az'])

        # Time shift to shift data to reference time
        tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']
        seistoplot.times = [x + tshift - seis[0].stats.traveltimes[phase] for x in seistoplot.times()]

        # Retrieve, filter and differentiate synthetics if wanted
        if syn:
            try:
                seissyn = seis.select(channel='BXT')[0]
                seissyn.filter('highpass', freq=fmin, corners=2, zerophase=True)
                seissyn.filter('highpass', freq=fmaz, corners=2, zerophase=True)
                seissyn.data = np.gradient(seissyn.data, seissyn.stats.delta)
                seissyn.data = seissyn.data / norm + np.round(seis[0].stats['az'])
                seissyn.times = [x - seis[0].stats.traveltimes[phase] for x in seissyn.times()]
            except:
                print('No synthetics available')

        # Plot data and synthetics if required
        plt.plot(seistoplot.times, seistoplot.data, 'k', linewidth=1)
        if syn:
            try:
                plt.plot(seissyn.times, seissyn.data, 'b', linewidth=1)
            except:
                print('No synthetics available')

        # Plot travel time predictions
        for k in seis[0].stats.traveltimes.keys():
            if seis[0].stats.traveltimes[k] != None:
                plt.plot(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az']), 'g', marker='o', markersize=4)
                # Uncomment line below if wish to display phase names ***TO IMPLEMENT*** Question at start to question if wanted
                #plt.text(seis[0].stats.traveltimes[k] - seis[0].stats.traveltimes[phase], np.round(seis[0].stats['az'])-0.5, k, fontsize=6)

        if syn:
            try:
                del seissyn
            except:
                print('No synthetics available')
        del seis    
                    
    # Put labels on graphs
    initDist = 85
    for i in range(4):   
        plt.subplot(1,4,i+1)
        plt.title('S dist < %d' % (initDist + (5 * i)), fontsize=10)
        if (i == 0):
            plt.ylabel('Azimuth (dg)')
        plt.ylim(round(min(azis) - 1), round(max(azis) + 1))        
        plt.xlim([-20, 60])           
        plt.xlabel('Time around predicted arrival (s)', fontsize=6)
        plt.tick_params(axis='both', which='major', labelsize=6)
        if switch_yaxis:
            plt.gca().invert_yaxis()      

    # Save file and show plot
    if syn == True:
        plt.savefig('Plots/' + name + '/' + 'azimuth.pdf')
    elif syn == False:
        plt.savefig('Plots/' + name + '/' + 'azimuth_noSyn.pdf')

    plt.show()
