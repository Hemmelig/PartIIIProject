#!/usr/bin/env python3

import obspy, os.path, glob, sys, scipy, scipy.signal, shutil
from obspy import read
import matplotlib.pyplot as plt
import numpy as np
from haversine import haversine, bearing
import subprocess
from subprocess import call
from obspy.taup.taup import getTravelTimes

event = sys.argv[1]

latmin = 48
latmax = 70
lonmin = -167
lonmax = -130

dir = 'Data/' + event
seislist = glob.glob(dir + '/*PICKLE')
seis = read(seislist[0], format='PICKLE')

seisT = seis.select(channel='BHT')[0]

evdp = str(seis[0].stats['dist'])

eventLat = seis[0].stats['evla']
eventLon = seis[0].stats['evlo']

outfile = 'receivers_forcsem_' + event + '_distaz.dat'
f = open(outfile, 'w')

nw = 'AA'

for i in range((lonmax - lonmin + 1)):
    slon = lonmin + (i * 1)
    if i < 10:
        ii = str(0) + str(i)
    else:
        ii = str(i)
    for j in range((latmax - latmin + 1)):
        slat = latmin + (j * 1)
        if j < 10:
            jj = str(0) + str(j)
        else:
            jj = str(j)

        sta = ii + jj

        dist = haversine(eventLat, eventLon, slat, slon)
        az   = bearing(eventLat, eventLon, slat, slon)

        f.write('%s %s %2.3f %3.3f\n' % (nw, sta, dist, az))

        synTc = seisT.copy()
        synTc.stats['channel'] = 'syncsemcomp'

        # Load synthetic
        filename = dir + 'CSEM/UT_' + nw + '_' + sta
        f2 = open(filename, 'r')
        timesyn = []
        sync = []
        lines = f2.readlines()
        for line in lines:
            line = line.split()
            timesyn.append(float(line[0]))
            sync.append(float(line[1]))

        # Setting the time array and the data, and stats it needs for filtering
        synTc.time = np.array(timesyn)
        synTc.data = np.array(sync)
        synTc.stats.npts = len(synTc.data)
        synTc.stats.delta = synTc.time[1] - synTc.time[0]
        synTc.stats.sampling_rate = 1. / synTc.stats.delta

        if not hasattr(synTc.stats, 'traveltimes'):
            synTc.stats.traveltimes = dict()

        # Add travel time predictions
        test = ['taup_time -mod prem -deg ' + str(dist) + ' -h ' + evdp + ' -ph S']
        out = subprocess.check_output(test, shell=True, universal_newlines=True)
        t = out.split()
        print(t)
        l = [x for x in range(len(t)) if t[x] == 'S']
        try:
            time = float(t[l[0] + 1])
            print('S', time)
        except:
            time = None
        synTc.stats.traveltimes['S'] = time

        f2.close()

        filename = dir + '/Synthetics/UT_' + nw + '_' + sta
        
        synTc.write(filename, format='PICKLE')

f.close()
