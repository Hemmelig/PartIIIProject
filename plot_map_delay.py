#!/usr/bin/env python3

import sys, obspy, obspy.signal, os.path, glob, scipy, subprocess
sys.path.append('/raid2/sc845/Python/lib/python/mpl_toolkits/')
sys.path.append('/raid2/sc845/Python/lib/python/')
sys.path.append('/raid2/sc845/Tomographic_models/LMClust/')
import LMClust_g6_ryb, mpl_toolkits, mpl_toolkits.basemap
from obspy import read
from obspy.core import Stream
import matplotlib.pyplot as plt
import numpy as np
print(mpl_toolkits.basemap.__path__)
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.image as mpimg
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

# Event to plot
name = sys.argv[1]

# Make directories for data
dir = 'Plots'
if not os.path.exists(dir):
    os.makedirs(dir)
dir = dir + '/' + name + '/'
if not os.path.exists(dir):
    os.makedirs(dir)

# Map parameters
plot_background = False
depth_background = 2700. # Clustering analysis is only down to 2700 km
midlat = 17.5
midlon = 192.5

dir = 'Data/' + name + '/'
seislist = glob.glob(dir + '/*PICKLE')

# Read in data from textfile
f = open(dir + 'peakData.txt', 'r')
lines = f.readlines()
delayTimes = []
amplRatios = []
azis = []
dists = []
for x in lines:
    delayTimes.append(x.split(' ')[0])
    amplRatios.append(x.split(' ')[1])
    azis.append(x.split(' ')[2])
    dists.append(x.split(' ')[3])    
f.close()

print(delayTimes, amplRatios, azis, dists)

# Start map
m = Basemap(projection='ortho', lat_0=midlat, lon_0=midlon, resolution='i')
clip_path = m.drawmapboundary()

# Plot background from the cluster analysis of Cottaar &Lekic 2016
if plot_background:
    LMC = LMClust_g6_ryb.LMClust()
    LMC.read('/raid2/sc845/Tomographic_models/LMClust/', 'clustgr.txt')
    lon, lat, layer = LMC.get_slice(depth_background)
    x, y = m(lon.ravel(), lat.ravel())
    x = x.reshape(np.shape(layer))
    y = y.reshape(np.shape(layer))
    minval = np.min(np.min(layer))
    maxval = np.max(np.max(layer)) + .1
    m.pcolor(x, y, layer, cmap=LMC.rgb_map_light, linewidth=0, rasterized=True)
    
m.drawcoastlines()

# Loop through stations
plot_event = True
for s in range(len(seislist)):
    print(s, seislist[s])
    seis = read(seislist[s], format='PICKLE')

    # Get event-station pair delay time and amplitude ratio
    delayTime = delayTimes[s]
    amplRatio = amplRatios[s]
    print(delayTime)

    elat = seis[0].stats['evla']
    elon = seis[0].stats['evlo']
    depth = seis[0].stats['evdp']

    # Plot event star for first station
    if plot_event:
        x2, y2 = m(elon+360., elat)
        m.scatter(x2, y2, s=265, marker='*', facecolors='y', alpha=1)
        plot_event = False

    # Plot station location
    slat = seis[0].stats['stla']
    slon = seis[0].stats['stlo']
    x2, y2 = m(slon, slat)
    m.scatter(x2, y2, s=35, c='w', marker='^', alpha=1)

    # Get bounce location of ScS phase
    turnloc = seis[0].stats.turnpoints["ScS"]
    tlon = turnloc[2]
    x2, y2 = m(tlon, turnloc[1])
    print(x2, y2)
    if (delayTime == 0):
        m.scatter(x2, y2, s=35, marker='x', alpha=1)
    else:
        m.scatter(x2, y2, s=35, c=delayTime, marker='o', alpha=1)

plt.colorbar()

# Location of ULVZ
x0, y0 = -167.5, 17.5
R = 7.5

m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.5)

plt.title('Turning depths in km')
plt.savefig('Plots/' + name + '/' + 'delayTimeMap.png')
plt.savefig('Plots/' + name + '/' + 'delayTimeMap.pdf')
            
plt.show()
