#!/usr/bin/env python3

import sys
sys.path.append('/raid2/sc845/Python/lib/python/mpl_toolkits/')
sys.path.append('/raid2/sc845/Python/lib/python/')
sys.path.append('/raid2/sc845/Tomographic_models/LMClust/')
import LMClust_g6_ryb
import obspy
from obspy import read
from obspy.core import Stream
import obspy.signal
import matplotlib.pyplot as plt
import os.path
import glob
import numpy as np
import scipy
import mpl_toolkits
import mpl_toolkits.basemap
print(mpl_toolkits.basemap.__path__)
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
import matplotlib.image as mpimg
import matplotlib.cm as cm
import subprocess
from matplotlib.colors import LinearSegmentedColormap

events = ["20170224", "20160528", "20160924"]

# Map parameters
plot_background = True
depth_background = 2700. # Clustering analysis is only down to 2700 km
midlat = 17.5
midlon = 192.5

fig = plt.figure(figsize=(6,10))

# Start map
m = Basemap(projection='ortho', lat_0=midlat, lon_0=midlon, resolution='f', llcrnrx=-2500000, llcrnry=-5000000, urcrnrx=3500000, urcrnry=5500000)

m.drawmapboundary()
cm = plt.cm.get_cmap('jet')

# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))

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

for i in range(len(events)):
    plot_event = True
    event = events[i]

    dir = 'Data/' + event + '/'

    f = open(dir + 'piercePoints.txt', 'r')
    lines = f.readlines()

    slats = []
    slons = []
    plat1 = []
    plon1 = []
    plat2 = []
    plon2 = []

    for x in lines:
        tmp = x.split(':')
        if plot_event:
            elat = float(tmp[0])
            elon = float(tmp[1])
            plot_event = False
        slats.append(float(tmp[2]))
        slons.append(float(tmp[3]))
        plat1.append(float(tmp[4]))
        plon1.append(float(tmp[5]))
        plat2.append(float(tmp[6]))
        plon2.append(float(tmp[7]))

    for s in range(len(slats)):
        print(s + 1, len(slats))

        m.drawgreatcircle(elon, elat, slons[s], slats[s], linewidth=0.5, color='gray', alpha=0.5, zorder=1)
        m.drawgreatcircle(plon1[s], plat1[s], plon2[s], plat2[s], linewidth=0.5, color='g', zorder=2)
        x1, y1 = m(slons[s], slats[s])
        m.scatter(x1, y1, s=40, marker='^', facecolors='white', alpha=1, zorder=2)

    x, y = m(elon + 360., elat)

    if (i == 2):
        m.scatter(x, y, s=265, marker='*', facecolors='y', alpha=1, zorder=2, label='Event')
        m.scatter(x1, y1, s=40, marker='^', facecolors='white', alpha=1, zorder=2, label='Station')
    else:
        m.scatter(x, y, s=265, marker='*', facecolors='y', alpha=1, zorder=2)


x0, y0 = -167.5, 17.5
R = 7.5

m.tissot(x0, y0, R, 100, facecolor='b', alpha=0.4, label='ULVZ')

plt.legend(loc=4)
         
fig.savefig('Plots/eventMap.png', bbox_inches='tight')
fig.savefig('Plots/eventMap.pdf', bbox_inches='tight')
         
plt.show()
