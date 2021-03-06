#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy, obspy.signal, os.path, time, glob, shutil, scipy, subprocess, sys
from obspy import read, UTCDateTime
from obspy.core import Stream, event, trace
from obspy.taup.taup import getTravelTimes
import matplotlib.pyplot as plt
import numpy as np
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call
import iterativeDecon as itd

#####################################################################################

def deconvolve_data(dataType, name, sourceName, idx1, idx2, deconType):

    fmax = .1
    fmin = 0.05
    savedist = []
    
    dir = 'Data/' + str(name) + '/'
    if (dataType == 1):
        dir = dir + 'RealData/'
    if (dataType == 2):
        dir = dir + 'CSEM/'
    if (dataType == 3):
        dir = dir + 'AxiSEM/'

    sourceseis = read(dir + sourceName + '.PICKLE', format='PICKLE')
    Stime = sourceseis[0].stats.traveltimes['S']
    source = sourceseis.select(channel='BHT')[0]
    source.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    source.data = np.gradient(source.data, source.stats.delta)

    tshift = sourceseis[0].stats['starttime'] - sourceseis[0].stats['eventtime']
    source.time = source.times()
    source.time = source.time[idx1:idx2]
    source.data = source.data[idx1:idx2]    
    
    norm = 0.3 * np.max([np.max(source.data), np.max(source.data)])
    source.taper(max_percentage=0.05, type='cosine')

    seislist = glob.glob(dir + '/*PICKLE')

    plt.figure(figsize=(14,10))

    plt.subplot(1,3,1)
    plt.plot(source.time, source.data, norm, 'k')

    plt.xlim([-25., 150])
    plt.ylim([-(1 / 0.3), (1 / 0.3)])
    plt.title('Source waveform')
    plt.xlabel('Time around S arrival (s)')
    
    dirLength = len(dir)    
    
    savedir = dir + 'Deconvolved/'


    for i in range(len(seislist)):
        print(i + 1, len(seislist))

        station = seislist[i]
        station = station[dirLength:]
        station = station[:len(station) - 7]
        
        seis = read(seislist[i], format='PICKLE')
        Stime = seis[0].stats.traveltimes['S']
        SHref = seis.select(channel='BHT')[0]
        dist = seis[0].stats['dist']

        SHref.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        SHref.data = np.gradient(SHref.data, SHref.stats.delta)

        norm = 0.3 * np.max(abs(SHref.data))

        # Filter seismograms
        SHref = SHref.slice(seis[0].stats['eventtime'] + Stime - 25., seis[0].stats['eventtime'] + Stime + 150)
        RF = trace.Trace()
        
        if (deconType == 1):
            #RF.data, fit = itd.water_level_decon(SHref.data, source.data, 3.e-2, source.stats['delta'], 'cosine', .5, 25.)
            pass
        if (deconType == 2):
            RF.data, fit = itd.iterative_deconvolution(SHref.data, source.data, 200, source.stats['delta'], 'cosine', .5, 25.)
        
        RF.data = 3. * RF.data / np.max(RF.data)
        RF.times = SHref.times()
        print(RF.data)
        print(fit)
        dist = np.round(dist)
        plt.subplot(1,3,2)
        plt.plot(SHref.times() - 25., RF.data + dist, 'k')
        plt.fill_between(SHref.times() - 25., dist, RF.data + dist, where = RF.data + dist > dist, facecolor = 'r')
        plt.fill_between(SHref.times() - 25., dist, RF.data + dist, where = dist > RF.data + dist, facecolor = 'b')
        plt.subplot(1,3,3)
        plt.plot(SHref.times() - 25., SHref.data / norm + dist, 'k')
        plt.fill_between(SHref.times() - 25., dist, SHref.data / norm + dist, where = SHref.data / norm + dist > dist, facecolor = 'r')
        plt.fill_between(SHref.times() - 25., dist, SHref.data / norm + dist, where = dist > SHref.data / norm + dist, facecolor = 'b')
        savedist.append(dist)

        # Save deconvolved waveforms so they can be accessed again.
        saveName = savedir + station
        RF.write(saveName + '.PICKLE', format='PICKLE')
                


    plt.subplot(1,3,2)
    plt.xlim([-25, 150])
    plt.ylim([min(savedist) - 2, max(savedist) + 2])
    plt.title('Deconvolved waveforms')
    plt.xlabel('Time around S arrival (s)')
    plt.ylabel('Distance (dg)')
    plt.subplot(1,3,3)
    plt.xlim([-25, 150])
    plt.ylim([min(savedist) - 2, max(savedist) + 2])
    plt.xlabel('Time around S arrival (s)')
    plt.ylabel('Distance (dg)')

    plt.title('Waveforms ' + savedir)

    plt.savefig('Plots/' + str(name) + '/' + sourceName + '_deconvolved.png', bbox_inches='tight')
    plt.savefig('Plots/' + str(name) + '/' + sourceName + '_deconvolved.pdf', bbox_inches='tight')
    plt.show()
