#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy, obspy.signal, os.path, time, glob, shutil, scipy, subprocess, sys
from obspy import read, UTCDateTime
from obspy.core import Stream, event
from obspy.taup.taup import getTravelTimes
import matplotlib.pyplot as plt
import numpy as np
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call

#####################################################################################

def addEventCMT(name):
    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    with open(dir + '/cmtsource.txt', 'r') as inf:
        srcdict = eval(inf.read())

    for s in range(len(seislist)):
        print(s)

        seis = read(seislist[s], format='PICKLE')

        # Event time
        time = UTCDateTime(srcdict['year'], srcdict['month'], srcdict['day'], srcdict['hour'], srcdict['min'], srcdict['sec'], srcdict['msec'])
        seis[0].stats['eventtime'] = time
        seis[1].stats['eventtime'] = time
        seis[2].stats['eventtime'] = time

        # Event location
        seis[0].stats['evla'] = srcdict['latitude']
        seis[0].stats['evlo'] = srcdict['longitude']

        # Event depth
        seis[0].stats['evdp'] = srcdict['depth']

        seis.write(seislist[s], format='PICKLE')

    print('CMT data for this event has been successfully added')
