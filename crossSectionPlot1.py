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
from scipy.interpolate import spline
import matplotlib.gridspec as gridspec

# Draw transect markers

def draw_screen_poly(lats, lons, m):
    x, y = m(lons, lats)
    xy = list(zip(x, y))
    poly = Polygon(xy, facecolor='red', alpha=0.2)
    plt.gca().add_patch(poly)

# Event to plot
names = ['20170224', '20160528']
pCorrections = [7.3, 7.6]

# Initialise figure
fig = plt.figure(figsize=(12,3))
gs = gridspec.GridSpec(1, 3)
ax = fig.add_subplot(gs[0], )

# Draw basemap
m = Basemap(llcrnrlon=-176, llcrnrlat=10, urcrnrlon=-157, urcrnrlat=26, resolution='f')
m.drawmapboundary(fill_color='#7777ff')
m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
m.drawcoastlines()

cm = plt.cm.get_cmap('jet')

# Location of ULVZ
x0, y0 = -167.5, 17.5
R = 7.5
m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.2)

lats1 = [17, 23, 23, 17]
lons1 = [-170, -170, -166, -166]

lats2 = [17.5, 21.5, 21.5, 17.5]
lons2 = [-172, -172, -166, -166]

draw_screen_poly(lats1, lons1, m)
draw_screen_poly(lats2, lons2, m)

# draw parallels.
parallels = np.arange(0.,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=16)
# draw meridians
meridians = np.arange(180.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=16)

tlons = []
tlats = []

for i in range(len(names)):
    name = names[i]
    pCorrection = pCorrections[i]
    dir = 'Data/' + name + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    # Read in data from textfile for delay times
    f = open(dir + 'peakData.txt', 'r')
    lines = f.readlines()
    delayTimes = []
    xs = []
    ys = []
    xs2 = []
    ys2 = []
    dtdt = []

    for x in lines:
        delayTimes.append(float(x.split(':')[0]))
    f.close()

    # Loop through stations
    for s in range(len(seislist)):
        print(s + 1, len(seislist), seislist[s])
        seis = read(seislist[s], format='PICKLE')

        # Get event-station pair delay time and amplitude ratio
        delayTime = delayTimes[s]

        # Get bounce location of ScS phase
        turnloc = seis[0].stats.turnpoints["ScS"]
        tlon = turnloc[2]
        tlat = turnloc[1]
        x, y = m(tlon, tlat)

        if (i == 0):
            dtPred = seis[0].stats.traveltimes["ScS"] - seis[0].stats.traveltimes["S"]
        else:
            dtPred = 0.
    
        if (delayTime == 0):
            xs2.append(x)
            ys2.append(y)
        else:
            xs.append(x)
            ys.append(y)
            tlons.append(tlon)
            tlats.append(tlat)
            dtdt.append(delayTime - dtPred - pCorrection)
            
    plot = m.scatter(xs, ys, s=50, c=dtdt, vmin=1, vmax=8, marker='o', alpha=1, cmap=cm, zorder=2)
    m.scatter(xs2, ys2, s=25, c='gray', marker='^', alpha=1, zorder=1)

# Section for producing heights plots
# Need to implement loop over 3 different shear wave velocity reductions

dvs = [24, 20]
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

dir1 = 'Data/'

for t in range(len(dvs)):
    dv = dvs[t]

    # Read in data from textfile for heights
    f = open(dir1 + 'convHeights_' + str(dv) + '.txt', 'r')
    lines = f.readlines()
    heights = []

    for x in lines:
        heights.append(float(x.split(':')[0]))
    f.close()    

    Lon_height = []
    Lon_height_std = []
    lonBin = [-172, -171.5, -171, -170.5, -170, -169.5, -169, -168.5, -168, -167.5, -167, -166.5, -166]
    for i in range(len(lonBin)):
        tmp_height = []
        for s in range(len(tlons)):
            if (tlons[s] >= lonBin[i] - 0.5 and tlons[s] <= lonBin[i] + 0.5 and tlats[s] >= 17.5 and tlats[s] <= 21.5):
                if (heights[s] == 0):
                    pass
                else:
                    tmp_height.append(heights[s])
        if not tmp_height:
            Lon_height.append(0)
            Lon_height_std.append(0)
        else:
            Lon_height.append(np.mean(tmp_height))
            Lon_height_std.append(np.std(tmp_height))
        del tmp_height

    indices = [i for i, x in enumerate(Lon_height) if x == 0]
    for i in range(len(indices)):
        j = indices[i]
        del Lon_height[j]
        del Lon_height_std[j]
        del lonBin[j]
        indices = [x - 1 for x in indices]

    # Longitude plot

    Lon_height1 = []
    Lon_height2 = []

    for i in range(len(Lon_height)):
        Lon_height1.append(Lon_height[i] + Lon_height_std[i])
        Lon_height2.append(Lon_height[i] - Lon_height_std[i])

    x_sm = np.array(lonBin)
    y_sm = np.array(Lon_height)
    y_sm1 = np.array(Lon_height1)
    y_sm2 = np.array(Lon_height2)

    x_smooth = np.linspace(x_sm.min(), x_sm.max(), 30)
    y_smooth = spline(lonBin, Lon_height, x_smooth)
    y_smooth1 = spline(lonBin, Lon_height1, x_smooth)
    y_smooth2 = spline(lonBin, Lon_height2, x_smooth)

    if (t == 0):
        ax1.plot(x_smooth, y_smooth, color='b', linewidth=1, label='-' + str(dv) + '%')
        ax1.plot(x_smooth, y_smooth1, 'b--',  linewidth=1)
        ax1.plot(x_smooth, y_smooth2, 'b--', linewidth=1)

        #ax1.fill_between(x_smooth, y_smooth2, y_smooth1, facecolor='blue', alpha=0.2)
    if (t == 1):
        ax1.plot(x_smooth, y_smooth, color='g', linewidth=1, label='-' + str(dv) + '%')
        ax1.plot(x_smooth, y_smooth1, 'g--',  linewidth=1)
        ax1.plot(x_smooth, y_smooth2, 'g--', linewidth=1)

        #ax1.fill_between(x_smooth, y_smooth2, y_smooth1, facecolor='green', alpha=0.2)        

    # Latitude plot
        
    Lat_height = []
    Lat_height_std = []
    latBin = [17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23]
    for i in range(len(latBin)):
        tmp_height = []
        for s in range(len(tlats)):
            if (tlats[s] >= latBin[i] - 0.5 and tlats[s] <= latBin[i] + 0.5 and tlons[s] >= -170 and tlons[s] <= -166):
                if (heights[s] == 0):
                    pass
                else:
                    tmp_height.append(heights[s])
        if not tmp_height:
            Lat_height.append(0)
            Lat_height_std.append(0)
        else:
            Lat_height.append(np.mean(tmp_height))
            Lat_height_std.append(np.std(tmp_height))
        del tmp_height

    indices = [i for i, x in enumerate(Lat_height) if x == 0]
    for i in range(len(indices)):
        j = indices[i]
        del Lat_height[j]
        del Lat_height_std[j]
        del latBin[j]
        indices = [x - 1 for x in indices]

    Lat_height1 = []
    Lat_height2 = []

    for i in range(len(Lat_height)):
        Lat_height1.append(Lat_height[i] + Lat_height_std[i])
        Lat_height2.append(Lat_height[i] - Lat_height_std[i])

    x_sm = np.array(latBin)
    y_sm = np.array(Lat_height)
    y_sm1 = np.array(Lat_height1)
    y_sm2 = np.array(Lat_height2)

    x_smooth = np.linspace(x_sm.min(), x_sm.max(), 30)
    y_smooth = spline(latBin, Lat_height, x_smooth)
    y_smooth1 = spline(latBin, Lat_height1, x_smooth)
    y_smooth2 = spline(latBin, Lat_height2, x_smooth)

    if (t == 0):
        ax2.plot(x_smooth, y_smooth, color='b', linewidth=1, label='-' + str(dv) + '%')
        ax2.plot(x_smooth, y_smooth1, 'b--',  linewidth=1)
        ax2.plot(x_smooth, y_smooth2, 'b--', linewidth=1)

        #ax2.fill_between(x_smooth, y_smooth2, y_smooth1, facecolor='blue', alpha=0.2)
    if (t == 1):
        ax2.plot(x_smooth, y_smooth, color='g', linewidth=1, label='-' + str(dv) + '%')
        ax2.plot(x_smooth, y_smooth1, 'g--',  linewidth=1)
        ax2.plot(x_smooth, y_smooth2, 'g--', linewidth=1)

        #ax2.fill_between(x_smooth, y_smooth2, y_smooth1, facecolor='green', alpha=0.2)    

# Configure axes

ax.set_title('a', fontsize=26, x=0.05, weight='bold')
ax1.set_title('b', fontsize=26, x=0.05, weight='bold')
ax2.set_title('c', fontsize=26, x=0.05, weight='bold')

ax1.set_ylabel('Height above CMB (km)', fontsize=18)
ax2.set_ylabel('Height above CMB (km)', fontsize=18)
ax1.set_xlabel('$^{\circ}$ Longitude', fontsize=18)
ax2.set_xlabel('$^{\circ}$ Latitude', fontsize=18)
ax1.tick_params(labelsize=16)
ax2.tick_params(labelsize=16)

# Colorbar axes
#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
#plt.colorbar(plot, cax=cbaxes)

ax1.set_xlim([-172, -166])
ax1.set_ylim([15, 45])
ax2.set_xlim([17, 23])
ax2.set_ylim([15, 45])

handles1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(handles1, labels1, loc=2)
handles2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles2, labels2, loc=2)

plt.show()

fig.savefig('Plots/heightMap.png', bbox_inches='tight')
fig.savefig('Plots/heightMap.pdf', bbox_inches='tight')
