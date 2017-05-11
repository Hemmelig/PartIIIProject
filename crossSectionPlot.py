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
from matplotlib.patches import Polygon

def draw_screen_poly(lats, lons, m):
    x, y = m(lons, lats)
    print(x, y)
    xy = list(zip(x, y))
    print(xy)
    poly = Polygon(xy, facecolor='red', alpha=0.2)
    plt.gca().add_patch(poly)

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

tlons = []
tlats = []

pCorrection = 14.6 / 2.
for x in lines:
    delayTimes.append(float(x.split(':')[0]))
    amplRatios.append(float(x.split(':')[1]))
    azis.append(float(x.split(':')[2]))
    dists.append(float(x.split(':')[3]))
    comments.append(x.split(':')[4])
f.close()

fig = plt.figure()
ax = fig.add_subplot(2,2,1)

m = Basemap(llcrnrlon=-176, llcrnrlat=8, urcrnrlon=-156, urcrnrlat=28, resolution='h')
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
for s in range(len(seislist)):
    print(s + 1, len(seislist[s]), seislist[s])
    seis = read(seislist[s], format='PICKLE')

    # Get event-station pair delay time and amplitude ratio
    delayTime = delayTimes[s]
    amplRatio = amplRatios[s]

    # Get bounce location of ScS phase
    turnloc = seis[0].stats.turnpoints["ScS"]
    tlon = turnloc[2]
    tlat = turnloc[1]
    x, y = m(tlon, tlat)

    dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
    
    if (delayTime == 0):
        xs2.append(x)
        ys2.append(y)
    else:
        xs.append(x)
        ys.append(y)
        tlons.append(tlon)
        tlats.append(tlat)
        dtdt.append(delayTime - dtPred - pCorrection)

lats1 = [15, 23, 23, 15]
lons1 = [-170, -170, -166, -166]

lats2 = [17.5, 21.5, 21.5, 17.5]
lons2 = [-173, -173, -164, -164]

draw_screen_poly(lats1, lons1, m)
draw_screen_poly(lats2, lons2, m)
        
plot = m.scatter(xs, ys, s=25, c=dtdt, vmin=1, vmax=8, marker='o', alpha=1, cmap=cm)

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(180.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

ax1 = fig.add_subplot(2,2,2)
Lon_tlons = [np.round(x) for x in tlons]
Lon_dtdt = []
Lon_dtdt_std = []
lonBin = [-173, -172, -171, -170, -169, -168, -167, -166, -165, -164]
for i in range(len(lonBin)):
    tmp_dtdt = []
    for s in range(len(Lon_tlons)):
        if (Lon_tlons[s] == lonBin[i] and tlats[s] >= 17.5 and tlats[s] <= 21.5):
            tmp_dtdt.append(dtdt[s])
    Lon_dtdt.append(np.mean(tmp_dtdt))
    Lon_dtdt_std.append(np.std(tmp_dtdt))
    del tmp_dtdt

ax1.errorbar(lonBin, Lon_dtdt, yerr=Lon_dtdt_std)

ax2 = fig.add_subplot(2,2,3)
Lat_tlons = [np.round(x) for x in tlats]
Lat_dtdt = []
Lat_dtdt_std = []
latBin = [15, 16, 17, 18, 19, 20, 21, 22, 23]
for i in range(len(latBin)):
    tmp_dtdt = []
    for s in range(len(Lat_tlons)):
        if (Lat_tlons[s] == latBin[i] and tlons[s] >= -170 and tlons[s] <= -166):
            tmp_dtdt.append(dtdt[s])
    Lat_dtdt.append(np.mean(tmp_dtdt))
    Lat_dtdt_std.append(np.std(tmp_dtdt))
    del tmp_dtdt

ax2.errorbar(latBin, Lat_dtdt, yerr=Lat_dtdt_std)

# Colorbar axes
#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
#plt.colorbar(plot, cax=cbaxes)

#plt.savefig('Plots/' + name + '/' + 'heightMap.png', bbox_inches='tight')
#plt.savefig('Plots/' + name + '/' + 'heightMap.pdf', bbox_inches='tight')
            
plt.show()
