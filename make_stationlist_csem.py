#!/usr/bin/env python3

#import obspy, obspy.signal, os.path, time, glob, shutil, scipy, subprocess, sys
import obspy, os.path, glob, sys
from obspy import read
#from obspy import UTCDateTime
#from obspy.core import Stream, event
#from obspy.taup.taup import getTravelTimes
#import matplotlib.pyplot as plt
#import numpy as np
#from obspy.xseed import Parser
#from obspy.arclink import Client as ARCLINKClient
#from obspy.fdsn import Client as IRISClient
#from subprocess import call

event = sys.argv[1]

outfile = 'receivers_forcsem_' + event + '.dat'
f = open(outfile, 'w')

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE') 
print(seislist)
azis = []

f.write('Number of stations is:\n')
f.write('%d \n' % len(seislist))
f.write('     nw stn lat lon:\n')
for s in range(len(seislist)):
       seis = read(seislist[s], format='PICKLE')
       sta  = seis[0].stats.station
       nw   = seis[0].stats.network
       slat = seis[0].stats['stla']
       slon = seis[0].stats['stlo']
       #print len(nw), len(sta)

       while len(nw) < 2:
              nw = '_' + nw
       while len(sta) < 4:
              sta = '_' + sta
       if len(sta) > 4:
              sta = sta[:4]
       print(len(nw), len(sta), nw, sta)

       f.write('%s %s  %2.3f  %3.3f\n'% (nw, sta, slat, slon))

f.close()
