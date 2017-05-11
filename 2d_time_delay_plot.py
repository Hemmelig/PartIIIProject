import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import sys
import csv

name = sys.argv[1]

# Defines array of distances to compute traveltime at
distances = [83, 84, 85, 86, 87, 88, 89, 90, 91]

# Anomalous layer
extent = (-0.25, -0.15, 5, 35)

dir = 'Data/PeakData/' + name

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

fig = plt.figure(figsize=(9,8))
fig.suptitle('2017-02-24', fontsize=20, x=0.47)
fig.text(0.47, 0.04, '% Shear velocity reduction', ha='center', fontsize=16)
fig.text(0.06, 0.5, 'Height of ULVZ', va='center', rotation='vertical', fontsize=16)

for i in range(len(distances)):
    with open(dir + '_' + str(int(distances[i])) + '_2d_dt.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        dt_arr = [[float(e) for e in r] for r in reader]

    ax = fig.add_subplot(3,3,i+1)
    ax.set_title(str(distances[i]) + r'$^{\circ}$', fontsize=8)
    cax = ax.imshow(dt_arr, interpolation='nearest', extent=extent, aspect='auto', vmin=1, vmax=8)

    ax.tick_params(labelsize=8)    
    
    if (i == 0 or i == 1 or i == 2 or i == 3 or i == 4 or i == 5):
        ax.set_xticklabels([])
    if (i == 1 or i == 2 or i == 4 or i == 5 or i == 7 or i == 8):
        ax.set_yticklabels([])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(cax, cax=cbar_ax)
cb.set_ticklabels([2, 3, 4, 5, 6, 7, 8])
cb.set_label('ScS delay time (s)', y=0.5)

plt.show()
        
fig.savefig('Plots/' + name + '/2d_dt_plots.png', bbox_inches='tight')
fig.savefig('Plots/' + name + '/2d_dt_plots.pdf', bbox_inches='tight')
