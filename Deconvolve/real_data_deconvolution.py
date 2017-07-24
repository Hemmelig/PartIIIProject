#!/usr/bin/env python3

import sys
import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import trace
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
from obspy import UTCDateTime
import receiver_function as rf
import subprocess

sourceDir = sys.argv[1]
dataDir   = sys.argv[2]

dir1 = 'Data/' + dataDir + '/'
seislist = glob.glob(dir1 + '*PICKLE')

dir2 = 'Data/' + sourceDir + '/'
sourcelist = glob.glob(dir2 + '*PICKLE')

c = 0

fmax = .5
fmin = 0.05 
savedist = []

sourceseis = read(sourcelist[1], format='PICKLE')
Ptime = sourceseis[0].stats.traveltimes['S']
SHref = sourceseis.select(channel='BHT')[0]
SHref.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
SHref.data = np.gradient(SHref.data, SHref.stats.delta)

tshift = sourceseis[0].stats['starttime'] - sourceseis[0].stats['eventtime']
print(sourceseis[0].stats['eventtime'], sourceseis[0].stats['starttime'])

source = SHref.slice(sourceseis[0].stats['eventtime'] + Ptime - 10., sourceseis[0].stats['eventtime'] + Ptime + 10.)
print(SHref.times(), source)
norm = 0.3 * np.max([np.max(source.data), np.max(source.data)])
source.taper(max_percentage=0.05, type='cosine')

plt.figure(figsize=(14,10))

plt.subplot(1,3,1)
plt.plot(SHref.times(), SHref.data / norm, 'k')

plt.xlim([Ptime - tshift - 10, Ptime - tshift + 10])
plt.ylim([-(1 / 0.3), (1 / 0.3)])
plt.title('Source waveform')
plt.xlabel('Time around S arrival (s)')


#for i in range(len(seislist)):
for i in range(10):
    print(i + 1, len(seislist))
            
    seis = read(seislist[i], format='PICKLE')
    Ptime = seis[0].stats.traveltimes['S']
    SHref = seis.select(channel='BHT')[0]
    dist = seis[0].stats['dist']
    c += 1
    print(Ptime)
    SHref.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    SHref.data = np.gradient(SHref.data, SHref.stats.delta)
            
    norm = 0.3 * np.max(abs(SHref.data))

    # Filter seismograms
    SHref = SHref.slice(seis[0].stats['eventtime'] + Ptime - 25., seis[0].stats['eventtime'] + Ptime + 150)    
    RF = trace.Trace()
    RF.data, fit = rf.water_level_decon(SHref.data, source.data, 3.e-2, source.stats['delta'], 'cosine', .5, 25.)
    RF.data = 3. * RF.data / np.max(RF.data)
    print(RF.data)
    print(fit)
    dist = np.round(dist)
    plt.subplot(1,3,2)
    plt.plot(SHref.times() - 25., RF.data + dist, 'k')
    plt.fill_between(SHref.times() - 25., dist, RF.data + dist, where = RF.data + dist > dist, facecolor = 'r')
    plt.fill_between(SHref.times() - 25., dist, RF.data + dist, where = dist > RF.data + dist, facecolor = 'b')         
    plt.subplot(1,3,3)
    plt.plot(SHref.times() - 25., SHref.data / norm + dist, 'k')
    #plt.plot(Pref.times() - Ptime, Pref.data / norm + dist, 'r')
    plt.fill_between(SHref.times() - 25., dist, SHref.data / norm + dist, where = SHref.data / norm + dist > dist, facecolor = 'r')
    plt.fill_between(SHref.times() - 25., dist, SHref.data / norm + dist, where = dist > SHref.data / norm + dist, facecolor = 'b')
    savedist.append(dist)
        
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

plt.title('Waveforms ' + dataDir )

plt.savefig('deconvolvedSV_' + dataDir + '.pdf')
plt.show()
