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

# Event to plot
name = sys.argv[1]

# Phases to plot turning points
phase = []
for i in range(2,len(sys.argv)):
    phase.append(sys.argv[i])

xs = []
ys = []
turns = []

# Map parameters
plot_background = True
depth_background = 2700. # Clustering analysis is only down to 2700 km
midlat = 17.5
midlon = 192.5

dir = 'Data/' + name + '/'
seislist = glob.glob(dir + '/*PICKLE')

# Start map
m = Basemap(projection='ortho',lat_0=midlat,lon_0=midlon,resolution='i')
#m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax#,
#                resolution='i',projection='merc',lon_0=27.,lat_0=46.7)
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
   print(s + 1, len(seislist[s]), seislist[s])
   seis = read(seislist[s], format='PICKLE')
   
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
   slon = [slon+360. if slon<0 else slon]
   x2, y2 = m(slon, slat)
   m.scatter(x2, y2, s=35, c='w', marker='^', alpha=1)

   # Plot phase turning point (in color of depth)
   # Upgoing Sdiff phases is plotted in black
   for ph in range(len(phase)):
      try:
          turnloc = seis[0].stats.turnpoints[phase[ph]]
          print(phase[ph], turnloc)
          if phase[ph] == 'Sdiff':
              lon = [turnloc[2]+360 if turnloc[2] < 0 else turnloc[2]]
              x2, y2 = m(lon, turnloc[1])
              m.scatter(x2, y2, s=35, c='r', marker='o', alpha=1)
              
              lon = [turnloc[4]+360 if turnloc[4] < 0 else turnloc[4]]
              x2, y2 = m(lon, turnloc[3])
              m.scatter(x2,y2,s=35,c='k',marker='o',alpha=1)

          else:
              lon = [turnloc[2]+360 if turnloc[2] < 0 else turnloc[2]]
              x2, y2 = m(lon, turnloc[1])
              xs.append(x2)
              ys.append(y2)
              turns.append(turnloc[0])
      except:
          print('No ' + phase[ph] + ' at this station')

plot = m.scatter(xs, ys, s=35, c=turns, vmin=2400, vmax=2900, marker='o', alpha=1)
          
plt.colorbar()

x0, y0 = -167.5, 17.5
R = 7.5

m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.05)

plt.title('Turning depths in km')
         
plt.savefig('Plots/' + name + '/' + 'turnPointMap.png')
plt.savefig('Plots/' + name + '/' + 'turnPointMap.pdf')
         
plt.show()
