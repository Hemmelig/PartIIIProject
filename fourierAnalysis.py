#!/usr/bin/env python3

import obspy
import obspy.signal
import os.path
import glob
import scipy
import sys
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

event = sys.argv[1]

fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE')

seis = read(seislist[0], format='PICKLE')
seistoplot = seis.select(channel='BHT')[0]
seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)

seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

Fs = 150.
Ts = 1. / Fs
t = np.arange(0, 1, Ts)

n = len(seistoplot.data)
k = np.arange(n)
T = n / Fs
frq = k / T
#frq = frq[range(n / 2)]

Y = np.fft.fft(seistoplot.data)
#Y = Y[range(n / 2)]

fig, ax = plt.subplots(2, 1)
ax[0].plot(seistoplot.times(), seistoplot.data)
ax[0].set_xlim(950, 1050)
ax[1].plot(frq, abs(Y))
ax[1].set_xlim(148,152)

plt.show()
