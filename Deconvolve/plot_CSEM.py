#!/usr/bin/env python3

import obspy, os.path, glob, sys, scipy, scipy.signal
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

event = sys.argv[1]

latmin = 48
latmax = 70
lonmin = -167
lonmax = -130

londiff = lonmax - lonmin
latdiff = latmax - latmin

nw = 'AA'
stas  = []
dists = []
azis  = []

norm = None

infile = 'receivers_forcsem_' + event + '_distaz.dat'
f = open(infile, 'r')

for columns in (row.strip().split() for row in f):
    stas.append(str(columns[1]))
    dists.append(float(columns[2]))
    azis.append(float(columns[3]))

dir = 'Data/' + event + 'CSEM/'

for s in range(len(dists)):
    print(s + 1, len(dists))
    sta  = stas[s]
    dist = dists[s]
    az   = azis[s]
    station_name = dir + 'UT_%s_%s' % (nw, sta)

    times = []
    seis  = []

    if dist < 80:
        plt.subplot(1,4,1)
    elif dist < 85:
        plt.subplot(1,4,2)
    elif dist < 90:
        plt.subplot(1,4,3)
    elif dist < 95:
        plt.subplot(1,4,4)

    crs = open(station_name, 'r')
    for columns in (row.strip().split() for row in crs):
        times.append(float(columns[0]))
        seis.append(float(columns[1]))

    if norm == None:
        norm = 1. * np.max(seis)
    seis = [x / norm for x in seis]
    plt.plot(times, seis + np.round(az), 'k', linewidth=1)
    
    
# Put labels on graphs
initDist = 80
for i in range(4):
    plt.ylim(round(min(azis) - 1), round(max(azis) + 1))
    plt.xlim([1500, 2000])
    plt.subplot(1,4,i+1)
    plt.title('S dist < %d' % (initDist + (5 * i)), fontsize=10)
    if (i == 1):
        plt.ylabel('Azimuth (dg)')
    plt.xlabel('Time around predicted arrival (s)', fontsize=6)
    plt.tick_params(axis='both', which='major', labelsize=6)

# Save file and show plot
plt.savefig(event + '_csem_azimuth.pdf')

plt.show()
