#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy
from obspy import read
from obspy.core import Stream
from obspy.core import event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import obspy.signal
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
import scipy
from obspy.io.xseed import Parser
from obspy.clients.arclink import Client as ARCLINKClient
from obspy.clients.fdsn import Client as IRISClient
from subprocess import call
import subprocess
from obspy.taup.taup import getTravelTimes
import sys

#####################################################################################

def turningPoints(*args):
    # Unpack args - they come in the form (name, [phases])
    name = args[0]
    phase = args[1]

    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')
    print(seislist)

    # Loop through data and read
    for s in range(len(seislist)):
        seis = read(seislist[s], format='PICKLE')

        # Start turning point dictionary
        if not hasattr(seis[0].stats, 'turnpoints'):
            seis[0].stats.turnpoints = dict()

        # Loop through phases and call TauP_time to get piercepoints
        for ph in range(len(phase)):

            elon = seis[0].stats['evlo']
            elat = seis[0].stats['evla']
            slat = seis[0].stats['stla']
            slon = seis[0].stats['stlo']
            depth = seis[0].stats['evdp']
            test = ['taup_pierce -mod prem -evt ' + str(elat) + ' ' + str(elon) + ' -sta ' + str(slat) + ' ' + str(slon) + ' -h ' + str(depth) + ' -ph ' + phase[ph] + ' -turn ']

            out = subprocess.check_output(test, shell=True, universal_newlines=True)
            t = out.split()
            if len(t) > 0:
                if phase[ph] == 'Sdiff':
                    turn2lat = float(t[-2])
                    turn2lon = float(t[-1])
                    turn1lat = float(t[-7])
                    turn1lon = float(t[-6])
                    turndepth = float(t[-4])

                    seis[0].stats.turnpoints[phase[ph]] = [turndepth, turn1lat, turn1lon, turn2lat, turn2lon]

                else:
                    turn1lat = float(t[-2])
                    turn1lon = float(t[-1])
                    turndepth = float(t[-4])

                    seis[0].stats.turnpoints[phase[ph]] = [turndepth, turn1lat, turn1lon]

            else:
                seis[0].stats.turnpoints[phase[ph]] = None

        print(seis[0].stats.turnpoints)

        # Write out seismogram again
        seis.write(seislist[s], format='PICKLE')

    print('Turning points added successfully')
