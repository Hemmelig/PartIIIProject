#!/usr/bin/env python3

import subprocess
import glob
import sys
import obspy
import os.path
from obspy import read
from obspy.core import Stream
import obspy.signal

events = ["20170224", "20160528", "20160924"]
piercedepth = 2851

for i in range(len(events)):
    event = events[i]
    dir = 'Data/' + event + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    slats = []
    slons = []

    pierce1lat = []
    pierce1lon = []
    pierce2lat = []
    pierce2lon = []
    
    
    get_event = True

    for s in range(len(seislist)):
        print(s + 1, len(seislist), seislist[s])

        seis = read(seislist[s], format='PICKLE')

        if get_event:
            elat = seis[0].stats['evla']
            elon = seis[0].stats['evlo']
            depth = seis[0].stats['evdp']
            get_event = False

        slat = seis[0].stats['stla']
        slon = seis[0].stats['stlo']
        slats.append(slat)
        slons.append(slon)

        test = ['taup_pierce -mod prem -evt ' + str(elat) + ' ' + str(elon) + ' -sta ' + str(slat) + ' ' + str(slon) + ' -h ' + str(depth)  +'  -ph ScS -pierce ' + str(piercedepth)]

        out = subprocess.check_output(test, shell=True, universal_newlines=True)

        print(out)

        t = out.split()
        print(t)

        if len(t) > 0:
            l = [x for x in range(len(t)) if (t[x] == str(2851.0))]
            print(l)
            pierce1lat.append(float(t[l[0] + 2]))
            pierce1lon.append(float(t[l[0] + 3]))
            pierce2lat.append(float(t[l[1] + 2]))
            pierce2lon.append(float(t[l[1] + 3]))
        else:
            pierce1lat.append(float(0.))
            pierce1lon.append(float(0.))
            pierce2lat.append(float(0.))
            pierce2lon.append(float(0.))


    piercePoints = open(dir + 'piercePoints.txt', 'w')
    for j in range(len(slats)):
        piercePoints.write("%s:%s:%s:%s:%s:%s:%s:%s \n" % (elat, elon, slats[j], slons[j], pierce1lat[j], pierce1lon[j], pierce2lat[j], pierce2lon[j]))
    piercePoints.close()
