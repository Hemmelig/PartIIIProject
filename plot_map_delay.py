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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

# Event to plot
name = sys.argv[1]

# Map parameters
plot_background = False
depth_background = 2700. # Clustering analysis is only down to 2700 km
midlat = 17.5
midlon = 192.5

dir = 'Data/' + name + '/'
seislist = glob.glob(dir + '/*PICKLE')
print(len(seislist))

def onpick(event):
    ind = event.ind[0]
    print(ind)
    print('Azimuth: ', azis[ind])
    print('Distance: ', dists[ind])
    print('Comment: ', comments[ind])
    print('Delay time: ', delayTimes[ind])

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
pCorrection = 29.13 / 4.
for x in lines:
    delayTimes.append(float(x.split(':')[0]))
    amplRatios.append(float(x.split(':')[1]))
    azis.append(float(x.split(':')[2]))
    dists.append(float(x.split(':')[3]))
    comments.append(x.split(':')[4])
f.close()

fig = plt.figure()
ax = fig.add_subplot(111)

plt.title('Delay time map')

# Start map
m = Basemap(projection='cyl', llcrnrlon=-200, llcrnrlat=-30, urcrnrlon=-110, urcrnrlat=75, resolution='h')
m.drawmapboundary(fill_color='#7777ff')
m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
m.drawcoastlines()

# Location of ULVZ
x0, y0 = -167.5, 17.5
R = 7.5

m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.2)

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

# Loop through stations
plot_event = True
for s in range(len(seislist)):
    print(s + 1, seislist[s])
    seis = read(seislist[s], format='PICKLE')

    # Get event-station pair delay time and amplitude ratio
    delayTime = delayTimes[s]
    amplRatio = amplRatios[s]

    # Plot event star for first station
    if plot_event:
        elat = seis[0].stats['evla']
        elon = seis[0].stats['evlo']
        depth = seis[0].stats['evdp']        
        x2, y2 = m(elon, elat)
        m.scatter(x2, y2, s=265, marker='*', facecolors='y', alpha=1)
        plot_event = False

    # Plot station location
    slat = seis[0].stats['stla']
    slon = seis[0].stats['stlo']
    x2, y2 = m(slon, slat)
    m.scatter(x2, y2, s=30, c='w', marker='^', alpha=1)

    # Get bounce location of ScS phase
    turnloc = seis[0].stats.turnpoints["ScS"]
    tlon = turnloc[2]
    x3, y3 = m(tlon, turnloc[1])

    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    print(dtPred, delayTime)

    if (delayTime == 0):
        xs.append(x3)
        ys.append(y3)
        dtdt.append(0)
        #m2.scatter(x3, y3, s=45, c='gray', marker='o', alpha=1)
    else:    
        xs.append(x3)
        ys.append(y3)
        dtdt.append(delayTime - dtPred - pCorrection)

axins = zoomed_inset_axes(ax, 2, loc=4)
axins.set_xlim(-176, -156)
axins.set_ylim(8, 28)

plt.xticks(visible=False)
plt.yticks(visible=False)

m2 = Basemap(llcrnrlon=-176, llcrnrlat=8, urcrnrlon=-156, urcrnrlat=28, ax=axins, resolution='h')
m2.drawmapboundary(fill_color='#7777ff')
m2.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
m2.drawcoastlines()

plot = m2.scatter(xs, ys, s=25, c=dtdt, vmin=1, vmax=8, marker='o', alpha=1, cmap=cm, picker=True)
fig.canvas.mpl_connect('pick_event', onpick)

m2.tissot(x0, y0, R, 100, facecolor='g', alpha=0.05)

mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")

# Colorbar axes
cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
plt.colorbar(plot, cax=cbaxes)

plt.savefig('Plots/' + name + '/' + 'delayTimeMap.png', bbox_inches='tight')
plt.savefig('Plots/' + name + '/' + 'delayTimeMap.pdf', bbox_inches='tight')
            
plt.show()
