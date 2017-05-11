#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap

from obspy.imaging.beachball import beach

fname = 'Data/events.txt'

def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))

lonmin = 150
latmin = -40
lonmax = 200
latmax = -10

names = []
elons = []
elats = []
depths = []
focmecs = []

# Start map
#m = Basemap(projection='ortho', lat_0=midlat, lon_0=midlon, resolution='i')
m = Basemap(projection='merc', lon_0=177.5, lat_0=-25.5, resolution='i', llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax)
cm = plt.cm.get_cmap('jet')

m.drawmapboundary(fill_color='aqua')
# fill continents, set lake color same as ocean color.
m.fillcontinents(color='coral',lake_color='aqua')

# draw parallels and meridians.
# label parallels on right and top
# meridians on bottom and left
parallels = np.arange(-40.,-10.,10.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(150,200.,10.)
m.drawmeridians(meridians,labels=[True,False,False,True])

# Read in file with cmt data in it
with open(fname) as f:
    content = f.readlines()
    for x in content:
        foctmp = []
        var = x.split(',')
        names.append(str(var[0]))
        elats.append(float(var[1]))
        if (float(var[2]) < 0):
            elon = float(var[2]) + 360.
        else:
            elon = float(var[2])
        elons.append(elon)
        depths.append(float(var[3]))
        foctmp.append(float(var[4]))
        foctmp.append(float(var[5]))       
        foctmp.append(float(var[6]))      
        focmecs.append(foctmp)

colors = get_colors(depths, plt.cm.jet, 400, 650)
cm = plt.cm.get_cmap('jet')
        
x, y = m(elons, elats)
nlons = []
nlats = []
for i in range(len(names)):
    if (i == 3):
        nlons.append(elons[i] - 9.5)
        nlats.append(elats[i] + 0.5)
    else:
        nlons.append(elons[i] + 1.5)
        nlats.append(elats[i] + 0.5)       
x2, y2 = m(nlons, nlats)
ax = plt.gca()
for s in range(len(names)):
    #m.scatter(x[s], y[s], s=265, marker='*', facecolors='y', alpha=1)    
    b = beach(focmecs[s], xy=(x[s], y[s]), facecolor=colors[s], width=150000, linewidth=1)
    b.set_zorder(10)
    ax.add_collection(b)

    plt.annotate(s=names[s], xy=(x[s], y[s]))

plot = m.scatter(x, y, s=35, c=depths, vmin=400, vmax=650, marker='o', alpha=0.01, cmap=cm)
    
plt.colorbar(plot)
    
plt.show()
    
