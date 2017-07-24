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

dir2 = 'Data/' + sourceDir + '/'
sourcelist = glob.glob(dir2 + '*PICKLE')

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

# First plot - SHref plot
plt.subplot(1,2,1)
plt.plot(SHref.times(), SHref.data / norm, 'k')
plt.xlim([Ptime - tshift - 10, Ptime - tshift + 10])
plt.ylim([-(1 / 0.3), (1 / 0.3)])

timeLength = len(source.data)
testTimes = [0] * timeLength
for i in range(timeLength):
    testTimes[i] = i * SHref.stats.delta

plt.subplot(1,2,2)
plt.plot(testTimes, source.data, 'k')

plt.show()
