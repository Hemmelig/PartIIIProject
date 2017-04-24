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
comments = []
xs = []
ys = []
dtdt = []
for x in lines:
    delayTimes.append(float(x.split(':')[0]))
    amplRatios.append(float(x.split(':')[1]))
    azis.append(float(x.split(':')[2]))
    dists.append(float(x.split(':')[3]))
    comments.append(x.split(':')[4])
f.close()

# Start map
m = Basemap(projection='ortho', lat_0=midlat, lon_0=midlon, resolution='i')
clip_path = m.drawmapboundary()
cm = plt.cm.get_cmap('jet')

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

    # Plot event star for first station
    if plot_event:
        elat = seis[0].stats['evla']
        elon = seis[0].stats['evlo']
        depth = seis[0].stats['evdp']        
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
    x3, y3 = m(tlon, turnloc[1])

    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    print(dtPred, delayTime)

    if (delayTime == 0):
        m.scatter(x3, y3, s=35, c='k', marker='o', alpha=1)
    else:    
        xs.append(x3)
        ys.append(y3)
        dtdt.append(delayTime - dtPred)

plot = m.scatter(xs, ys, s=35, c=dtdt, vmin=11, vmax=20, marker='o', alpha=1, cmap=cm)

plt.colorbar(plot)

print(dtdt, max(dtdt))

# Location of ULVZ
x0, y0 = -167.5, 17.5
R = 7.5

m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.05)

plt.title('Delay time map')
plt.savefig('Plots/' + name + '/' + 'delayTimeMap.png')
plt.savefig('Plots/' + name + '/' + 'delayTimeMap.pdf')
            
plt.show()
