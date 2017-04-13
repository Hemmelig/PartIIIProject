#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy, obspy.signal, os.path, time, glob, shutil, scipy, subprocess, sys
from obspy import read
from obspy.core import Stream, event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call

#####################################################################################

def travelTimes(*args):
    # Unpack args - they come in the form (name, [phases])
    name = args[0]
    phase = args[1]
    

    # This code assumes the data is in a Data directory and in PICKLE format. Change here if different. 
    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE') 
    print(seislist)

    # Loop through data and read
    for s in range(len(seislist)):
        seis = read(seislist[s], format='PICKLE')

        # Start travel time dictionary
        if not hasattr(seis[0].stats,'traveltimes'):
            seis[0].stats.traveltimes = dict()

        # Loop through phases and call TauP_time to get travel times
        for ph in range(len(phase)):
            test = ['taup_time -mod prem -deg ' + str(seis[0].stats['dist']) + ' -h ' + str(seis[0].stats['evdp']) + ' -ph ' + phase[ph]]
            out = subprocess.check_output(test, shell=True, universal_newlines=True)
            t = out.split()
            print(t)
            l = [x for x in range(len(t)) if t[x] == phase[ph]]
            try:
                time = float(t[l[0] + 1])
                print(phase[ph], time)
            except:
                time = None
            seis[0].stats.traveltimes[phase[ph]] = time

        # Write out seismogram again
        seis.write(seislist[s], format='PICKLE')

    print('Travel times have been successfully added')
