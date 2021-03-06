#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy, obspy.signal, obspy.signal.rotate, os.path, time, glob, shutil, scipy, sys
from obspy import read
from obspy.core import Stream, event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np

#####################################################################################

def processData(name):
    # Data is assumed to be in PICKLE format in Data/<eventname>/Originals/
    dir = 'Data/' + name + '/'
    stalist = glob.glob(dir + '/Originals/*PICKLE')

    for s in range(len(stalist)):
        try:
            onestation = read(stalist[s], format='PICKLE')
            
            # Cut components to the same length
            onestation.trim(starttime=onestation[0].stats.starttime + 300, endtime=onestation[0].stats.endtime - 300)
            
            # Find if component names are different
            seisZ = onestation.select(channel='BHZ')
            if len(seisZ) > 1:
                seisZ = seisZ.select(location='10')

            seisN = onestation.select(channel='BHN')
            if len(seisN) == 0:
                seisN = onestation.select(channel='BH2')
            if len(seisN) > 1:
                seisN = seisN.select(location='10')

            seisE = onestation.select(channel='BHE')
            if len(seisE) == 0:
                seisE = onestation.select(channel='BH1')
            if len(seisE) > 1:
                seisE = seisE.select(location='10')

            print(seisZ, seisN, seisE)
            seisZ[0].stats = onestation[0].stats

            # Rotate components from North and East to Radial and Transverse
            [seisRtmp, seisTtmp] = obspy.signal.rotate.rotate_ne_rt(seisN[0].data, seisE[0].data, seisZ[0].stats['baz'])
            seisR = seisN[0].copy()
            seisR.stats['channel'] = 'BHR'
            seisR.data = seisRtmp
            seisT = seisN[0].copy()
            seisT.stats['channel'] = 'BHT'
            seisT.data = seisTtmp

            # Copy values into stats for vertical component
            seisZ[0].stats = onestation[0].stats

            # Produce new stream with Vertical, Radial and Transverse
            seisnew = Stream()
            seisnew.append(seisZ[0])
            seisnew.append(seisR)
            seisnew.append(seisT)
            # Write out in PICKLE format
            filename = dir + seisZ[0].stats.network + '.' + seisZ[0].stats.station + '.PICKLE'
            seisnew.write(filename, 'PICKLE')
        except:
            print('FAILED for', stalist[s])

    print('Data process complete')
