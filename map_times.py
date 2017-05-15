#!/usr/bin/env python3

import sys, obspy, obspy.signal, os.path, glob, scipy, subprocess
from obspy import read
from obspy.core import Stream
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.image as mpimg
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

# Event to plot
events = ['20170224', '20160528', '20160924']
pCorrections = [7.3, 7.6, 9.25]
titles = ['a', 'b', 'c']

# Set up figure
fig = plt.figure(figsize=(6,3))

for i in range(len(events)):
    dir = 'Data/' + events[i] + '/'
    seislist = glob.glob(dir + '/*PICKLE')

    # Initialise arrays
    delayTimes = []
    xs = []
    ys = []
    xs2 = []
    ys2 = []
    dtdt = []
    tlons = []
    tlats = []

    pCorrection = pCorrections[i]
    
    # Read in data
    f = open(dir + 'peakData.txt', 'r')
    lines = f.readlines()
    for x in lines:
        delayTimes.append(float(x.split(':')[0]))
    f.close()

    # Create axes
    ax = fig.add_subplot(1,3,i+1)
    ax.set_title(titles[i], fontsize=24, x=0.05, weight='bold')

    # Draw basemap
    m = Basemap(llcrnrlon=-176, llcrnrlat=10, urcrnrlon=-157, urcrnrlat=26, resolution='f')
    m.drawmapboundary(fill_color='#7777ff')
    m.fillcontinents(color='#ddaa66', lake_color='#7777ff', zorder=0)
    m.drawcoastlines()

    # Add ULVZ
    x0, y0 = -167.5, 17.5
    R = 7.5
    m.tissot(x0, y0, R, 100, facecolor='g', alpha=0.2)

    cm = plt.cm.get_cmap('jet')

    # Loop through stations
    for s in range(len(seislist)):
        print(s + 1, len(seislist), seislist[s])
        seis = read(seislist[s], format='PICKLE')

        # Get event-staton pair delay time
        delayTime = delayTimes[s]

        # Get bounce location of ScS phase
        turnloc = seis[0].stats.turnpoints['ScS']
        tlon = turnloc[2]
        tlat = turnloc[1]
        x, y = m(tlon, tlat)

        if (i == 0):   
            dtPred = seis[0].stats.traveltimes['ScS'] - seis[0].stats.traveltimes['S']
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
    # draw parallels.
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
        
# Add details to plots
cbaxes = fig.add_axes([0.15, 0.25, 0.7, 0.03]) # setup colorbar axes.
cb = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
cb.set_label('Sc$^{\prime}$S - ScS delay time (s)', fontsize=24)
cb.ax.tick_params(labelsize=14)
            
plt.show()

fig.savefig('Plots/eventsdt.png', bbox_inches='tight')
fig.savefig('Plots/eventsdt.pdf', bbox_inches='tight')
