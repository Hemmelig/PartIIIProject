#!/usr/bin/env python3

import obspy
import os.path
import glob
import sys
from obspy import read

def makeStationList(name):
    dir = 'Data/StationLists/'    
    outfile = dir + 'receivers_forcsem_' + name + '.dat'
    f = open(outfile, 'w')

    dir = 'Data/' + name + '/'       
    seislist = glob.glob(dir + '/*PICKLE')
    print(seislist)
    azis = []

    f.write('Number of stations: \n')
    f.write('%d \n' % len(seislist))
    f.write('     nw stn lat lon:\n')
    for s in range(len(seislist)):
        seis = read(seislist[s], format='PICKLE')
        sta  = seis[0].stats.station
        nw   = seis[0].stats.network
        slat = seis[0].stats['stla']
        slon = seis[0].stats['stlo']

        while len(nw) < 2:
            nw = '_' + nw
        while len(sta) < 4:
            sta = '_' + sta
        if len(sta) > 4:
            sta = sta[:4]

        print(nw, sta)

        f.write('%s %s  %2.3f  %3.3f\n' % (nw, sta, slat, slon))

    f.close()
