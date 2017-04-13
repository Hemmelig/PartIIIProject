#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import instaseis, obspy, sys, glob, obspy.signal.rotate
from obspy import read
import matplotlib.pyplot as plt

#####################################################################################

def addSynthetics(name, clean):
    print(instaseis.__path__)

    # Input directory with PICKLE data as argument
    dir = 'Data/' + name + '/'

    # Load database with Greens functions
    db = instaseis.open_db("/raid2/sc845/Instaseis/DB/10s_PREM_ANI_FORCES/")

    # Directory needs to contain a CMT source text file
    with open(dir + 'cmtsource.txt', 'r') as inf:
        srcdict = eval(inf.read())

    # Unpack dictionary
    latitude = srcdict['latitude']
    longitude = srcdict['longitude']
    depth_in_m = srcdict['depth'] * 1.e3
    m_rr = srcdict['Mrr'] / 1E7
    m_tt = srcdict['Mtt']  / 1E7
    m_pp = srcdict['Mpp'] / 1E7
    m_rt = srcdict['Mrt'] / 1E7
    m_rp = srcdict['Mrp']/ 1E7
    m_tp = srcdict['Mtp'] / 1E7
    origin_time = obspy.UTCDateTime(srcdict['year'], srcdict['month'], srcdict['day'], srcdict['hour'], srcdict['min'], srcdict['sec'], srcdict['msec'])
    
    # Read in source
    source = instaseis.Source(latitude=latitude, longitude=longitude, depth_in_m=depth_in_m, m_rr=m_rr, m_tt=m_tt, m_pp=m_pp, m_rt=m_rt, m_rp=m_rp, m_tp=m_tp, origin_time=origin_time)

    # Read and loop through station list
    stalist = glob.glob(dir + '/*PICKLE')

    for s in range(len(stalist)):
        seis = read(stalist[s], format='PICKLE')

        for tr in seis.select(channel='BX*'):
            seis.remove(tr)

        # Fixing some old mistake (not sure if still present, check with Sanne)
        seis[0].stats['channel'] = 'BHZ'
        print(seis[0].stats['channel'])
        receiver = instaseis.Receiver(latitude=seis[0].stats['stla'], longitude=seis[0].stats['stlo'], network=seis[0].stats['network'], station=seis[0].stats['station'])
        
        start = seis[0].stats['starttime']
        end = seis[0].stats['endtime']
        st = db.get_seismograms(source=source, receiver=receiver, kind='displacement', dt=0.1)

        # Rotate synthetics
        stE = st.select(channel='BXE')
        stN = st.select(channel='BXN')
        stZ = st.select(channel='BXZ')

        [stRtmp, stTtmp] = obspy.signal.rotate.rotate_ne_rt(stN[0].data, stE[0].data, seis[0].stats['baz'])

        stR = stN[0].copy()
        stR.stats['channel'] = 'BXR'
        stR.data = stRtmp
        stT = stN[0].copy()
        stT.stats['channel'] = 'BXT'
        stT.data = stTtmp

        seis += stR
        seis += stT
        seis += stZ

        for x in seis:
            print(x.stats['channel'])

        # Overwrites previous PICKLE with synthetics included
        seis.write(stalist[s], format='PICKLE')

    print('Synthetics for this event have been successfully added')
