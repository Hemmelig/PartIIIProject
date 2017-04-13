#!/usr/bin/env python3

#####################################################################################
### Import modules
#####################################################################################

import obspy
import os.path
import time
import obspy.geodetics.base
from obspy.clients.fdsn import Client as IRISClient
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np

#####################################################################################

def downloadData(name, latitude, longitude, st, et, maxrad, minmag, maxmag, distmin, distmax, azmin, azmax, lengthoftrace):
    # Load IRIS client
    irisclient = IRISClient("IRIS")

    # Just to ensure it is functioning
    print(name, latitude, longitude, st, et, maxrad, minmag, maxmag, distmin, distmax, azmin, azmax, lengthoftrace)
    
    # Convert time arguments
    starttime = UTCDateTime(st)
    endtime = UTCDateTime(et)

    # Define a filter band to prevent amplifying noise during the deconvolution
    # Currently not in use, as responses are not being removed
    fl1 = 0.005
    fl2 = 0.01
    fl3 = .5
    fl4 = 2.

    # Make directories for data
    dir = 'Data'
    if not os.path.exists(dir):
        os.makedirs(dir)
    dir = dir + '/' + name
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Directory for response files
    dirresp = dir + '/Responsefiles'
    if not os.path.exists(dirresp):
        os.makedirs(dirresp)

    # Directory for originals
    dir = dir + '/Originals'
    if not os.path.exists(dir):
        os.makedirs(dir)

    # Directory for plots (later)
    dir2 = 'Plots'
    if not os.path.exists(dir2):
        os.makedirs(dir2)
    dir2 = dir2 + '/' + name
    if not os.path.exists(dir2):
        os.makedirs(dir2)
   
    # Get event catalogue via IRIS client and save event parameters
    cat = irisclient.get_events(latitude=latitude, longitude=longitude, maxradius=maxrad, starttime=starttime, endtime=endtime, minmagnitude=minmag)
    evtlatitude = cat[0].origins[0]['latitude']
    evtlongitude = cat[0].origins[0]['longitude']
    evtdepth = cat[0].origins[0]['depth'] / 1.e3 # Convert to km from m
    evstarttime = cat[0].origins[0].time
    evendtime = evstarttime + lengthoftrace

    # Check how many events have been collected - if more than one, terminate script
    if len(cat) > 1:
        print('more than one event is selected')
        sys.exit()

    # Select what stations are present at the time
    inventory = irisclient.get_stations(starttime=evstarttime, endtime=evendtime)# Channel information does not seem to work ...,channel='BH*')
    count = 0
    inventory.get_response

    # Loop through networks
    for nw in inventory:
        print(nw)  

        # Loop through each station for each network
        for sta in nw:
            print(sta)
    
            seis = []
            distm, az, baz = obspy.geodetics.base.gps2dist_azimuth(evtlatitude, evtlongitude, sta.latitude, sta.longitude)
            # Convert distance from units of metres to degrees
            distdg = distm / (6371.e3 * np.pi / 180.)
            print(distdg, az)
            
            # Select data depending on distance and azimuths
            if distdg > distmin and distdg < distmax and az > azmin and az < azmax:
                # Try to download data
                try:
                    print('Trying', nw.code, sta.code, "*", "BH*", evstarttime, evendtime)
                    seis = irisclient.get_waveforms(nw.code, sta.code, "*", "BH*", evstarttime, evendtime, attach_response=True)
                except:
                    print('Failed',nw.code, sta.code, "*", "BH*", evstarttime, evstarttime + lengthoftrace)

            # If data has been found, save with file format PICKLE        
            if len(seis) > 1:
                # Get response file
                try:
                    seis.remove_response(output='DISP', pre_filt=[fl1,fl2,fl3,fl4])
	        
                    seis[0].stats['evla'] = evtlatitude
                    seis[0].stats['evlo'] = evtlongitude
                    seis[0].stats['evdp'] = evtdepth
                    seis[0].stats['stla'] = sta.latitude
                    seis[0].stats['stlo'] = sta.longitude
	
                    seis[0].stats['dist'] = distdg
                    seis[0].stats['az'] = az
                    seis[0].stats['baz'] = baz
		
                    seis[0].stats['station'] = sta.code
                    seis[0].stats['network'] = nw.code
			
                    filename = dir + '/' + seis[0].stats.station + '.' + seis[0].stats.network + '.PICKLE'
                    print('Writing to ' , filename)
                    count = count + 1
                    seis.write(filename, format='PICKLE')
                except:
                    print('Response removal or saving failed for, ' + nw.code + ' ' + sta.code)
        
    print('Seismograms found for ' + str(count) + ' stations')
            
