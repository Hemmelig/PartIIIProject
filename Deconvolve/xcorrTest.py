#!/usr/bin/env python3

import obspy
import obspy.signal
import os.path
import glob
import scipy
import sys
from obspy import read
from obspy.signal.cross_correlation import xcorr
import matplotlib.pyplot as plt
import numpy as np

event = sys.argv[1]
window = int(sys.argv[2])

fmin = 0.033
fmax = 0.1
dir = '../Data/' + event + '/'
seislist = glob.glob(dir + '*PICKLE')

seis1 = read(seislist[0], format='PICKLE')
seis2 = read(seislist[2], format='PICKLE')

seistoplot1 = seis1.select(channel='BHT')[0]
seistoplot2 = seis2.select(channel='BHT')[0]

seistoplot1.filter('highpass', freq=fmin, corners=2, zerophase=True)
seistoplot1.filter('lowpass', freq=fmax, corners=2, zerophase=True)
seistoplot2.filter('highpass', freq=fmin, corners=2, zerophase=True)
seistoplot2.filter('lowpass', freq=fmax, corners=2, zerophase=True)
seistoplot1.data = np.gradient(seistoplot1.data, seistoplot1.stats.delta)
norm1 = np.max(abs(seistoplot1.data))
seistoplot1.data = seistoplot1.data / norm1

seistoplot2.data = np.gradient(seistoplot2.data, seistoplot2.stats.delta)
norm2 = np.max(abs(seistoplot2.data))
seistoplot2.data = seistoplot2.data / norm2
print(np.max(seistoplot1.data), np.max(seistoplot2.data))

shift, coe = xcorr(seistoplot1.data, seistoplot2.data, window)

print(len(seistoplot1), seistoplot1.stats.delta)

tshift1 = seis1[0].stats['starttime'] - seis1[0].stats['eventtime']
tshift2 = seis2[0].stats['starttime'] - seis2[0].stats['eventtime']

print(tshift1, tshift2)

print(shift, coe)
