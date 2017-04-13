#!/usr/bin/env python3

import obspy, os.path, time, glob, shutil, sys
from obspy import read, UTCDateTime
from obspy.core import Stream
import matplotlib.pyplot as plt
import numpy as np


'''
This script can be used to select data. 
It will loop through the seismograms and show the Radial component (red) and Transverse component (black). 
Close the window by clicking ctrl+w. 
Then choose if you want the seismogram by typing 'y' or 'n' and pressing enter. If you type anything else it will show the seismogram for again.
Discarded data will be moved to a directory called Dump (you can always move it back from there). 
'''


event = sys.argv[1]

#Also plot synthetics
syn = True

## Filter freq
fmin = 0.033
fmax = 0.1

dir = 'Data/' + event + '/'
seislist = glob.glob(dir + '/*PICKLE') 

# Create directory to dump data into
dirdump = dir + 'Dump'
if not os.path.exists(dirdump):
    os.makedirs(dirdump)

for s in range(len(seislist)):
    print(s, len(seislist), seislist[s])
    seis = read(seislist[s], format='PICKLE')

    answer = '' 
    while answer != 'y' and answer != 'n':

        #plt.ion()
        plt.subplot(211)
        
        if (answer != 'p'):
            seis.filter('highpass',freq=fmin,corners=2,zerophase=True)
            seis.filter('lowpass',freq=fmax,corners=2,zerophase=True)
            
        phase='Sdiff'
        if seis[0].stats.traveltimes[phase] is None:
            phase='S'
 
        tshift = seis[2].stats['starttime'] - seis[2].stats['eventtime']

        if (answer != 'p'):
            seisR = seis.select(channel='BHR')[0]
            seisT = seis.select(channel='BHT')[0]
            scaler = np.max(seisT.data)
            
        # plotting data. Note, both components are normalized by the maximum amplitude on the Transverse component
        plt.subplot(2,1,2)
        plt.plot(seisR.times() + tshift, seisR.data / scaler, 'r', linewidth=3)
        plt.subplot(2,1,1)       
        plt.plot(seisT.times() + tshift, seisT.data / scaler, 'k', linewidth=3)

        # Plot synthetic data
        if syn == True:
             synT = seis.select(channel='BXT')[0]
             if (answer != 'p'):
                 synT.filter('highpass', freq=fmin, corners=2, zerophase=True)
                 synT.filter('lowpass', freq=fmax, corners=2, zerophase=True)
             synR=seis.select(channel='BXR')[0]
             if (answer != 'p'):
                 synR.filter('highpass', freq=fmin, corners=2, zerophase=True)
                 synR.filter('lowpass', freq=fmax, corners=2, zerophase=True)
             plt.subplot(2,1,2)
             plt.plot(synR.times(), synR.data / scaler, color=[0.5,0.5,0.5])
             plt.subplot(2,1,1)
             plt.plot(synT.times(), synT.data / scaler, color=[0.5,0.5,0.5])
         
        # Plotting travel time predictions
        for k in seis[0].stats.traveltimes.keys():
               if seis[0].stats.traveltimes[k]!=None:
                   plt.plot(seis[0].stats.traveltimes[k],0.0,'g',marker='o',markersize=4)
                   plt.text(seis[0].stats.traveltimes[k],-0.5, k)
        plt.title(str(seis[0].stats['dist']) +'   ' + str(seis[0].stats['az']))
        plt.subplot(2,1,1)        
        plt.ylim([-1.0,1.0])
        plt.xlim([seis[0].stats.traveltimes[phase]-100,seis[0].stats.traveltimes[phase]+100])
        plt.subplot(2,1,2)        
        plt.ylim([-1.0,1.0])
        plt.xlim([seis[0].stats.traveltimes[phase]-100,seis[0].stats.traveltimes[phase]+100])        

        plt.show()

        answer = input('Want this seismogram (y/n/p)? ')
 
        if (answer == 'y'):
            print('data accepted')
        if (answer == 'n'):
            shutil.move(seislist[s], dirdump)
        if (answer == 'p'):
            seisR.data = seisR.data * -1.
            seisT.data = seisT.data * -1.

            seisZ = seis.select(channel='BHZ')
            if len(seisZ) > 1:
                seisZ = seisZ.select(location='10')
            seisZ[0].stats = seis[0].stats
            
            seisnew = Stream()
            seisnew.append(seisZ[0])
            seisnew.append(seisR)
            seisnew.append(seisT)
            # Write out in PICKLE format
            filename = seislist[s]
            seisnew.write(filename, 'PICKLE')
