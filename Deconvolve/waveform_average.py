#!/usr/bin/env python3

import obspy, obspy.signal, os.path, glob, scipy, sys
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

event = sys.argv[1]

def findAlign(timeArray, seisArray):
    iInd = np.searchsorted(timeArray, -20)
    fInd = np.searchsorted(timeArray, 60)
    windowSeis = []
    for i in range(iInd, fInd):
        windowSeis.append(seisArray[i])

    maxInd = np.argmax(windowSeis)
    return maxInd, iInd

def shift_left_once(lst):
    temp = lst[0]
    for index in range(len(lst) - 1):
        lst[index] = lst[index + 1]         
    lst[index + 1] = temp

def shift(lst, n):
    """Shifts the lst over by n indices

    >>> lst = [1, 2, 3, 4, 5]
    >>> shift_left(lst, 2)
    >>> lst
    [3, 4, 5, 1, 2]
    """
    assert (n >= 0), "n should be non-negative integer"
    for _ in range(n):
        shift_left_once(lst)

# Frequencies for filter
fmin = 0.033
fmax = 0.1

dir = '../Data/' + event + '/'
seislist = glob.glob(dir + '*PICKLE')

seisTimes   = []

fig = plt.figure(figsize=(4,12))

# Loop through seismograms
for s in range(len(seislist)):
    print(s + 1, len(seislist))
    seis = read(seislist[s], format='PICKLE')

    seistoplot = seis.select(channel='BHT')[0]

    seistoplot.filter('highpass', freq=fmin, corners=2, zerophase=True)
    seistoplot.filter('lowpass', freq=fmax, corners=2, zerophase=True)
    seistoplot.data = np.gradient(seistoplot.data, seistoplot.stats.delta)

    # Align peaks of S arrival on t = 0
    tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']
    seistoplot.times = [x + tshift - seis[0].stats.traveltimes['S'] for x in seistoplot.times()]
    alignInd, iInd = findAlign(seistoplot.times, seistoplot.data)
    print(alignInd, iInd)
    talign = seistoplot.times[alignInd + iInd]
    print(talign)
    seistoplot.times = [x - talign for x in seistoplot.times]

    if s == 0:
        seisTimes = seistoplot.times
        seisAverage = [0] * len(seistoplot.data)

    if alignInd < 0:
        print('Uh oh...')
    if alignInd > 0:
        shift(seistoplot.data, alignInd)
    if alignInd == 0:
        pass

    norm = np.max(np.abs(seistoplot.data))
    plt.plot(seistoplot.times, seistoplot.data / norm + (s * 2), 'k')
    
    #tmp = seisAverage
            
    #seisAverage = [sum(x) for x in zip(tmp, seistoplot.data)]

    #del tmp

#norm = np.max(np.abs(seisAverage))
    
#plt.plot(seistoplot.times, seisAverage / norm, 'k')

plt.show()


