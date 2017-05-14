import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import sys
import csv

name = sys.argv[1]

# Defines array of distances to compute traveltime at
distances = [80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97]

# Anomalous layer
extent = (-0.25, -0.15, 4, 45)

dir = 'Data/PeakData/' + name

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

fig = plt.figure(figsize=(9,8))
fig.text(0.47, 0.03, '% Shear velocity reduction', ha='center', fontsize=26)
fig.text(0.07, 0.5, 'Height of ULVZ (km)', va='center', rotation='vertical', fontsize=26)

for i in range(len(distances)):
    with open(dir + '_' + str(int(distances[i])) + '_2d_dt.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        dt_arr = [[float(e) for e in r] for r in reader]
        dt_arr = np.delete(dt_arr, 0, 1)

    ax = fig.add_subplot(3,6,i+1)
    ax.set_title(str(distances[i]) + r'$^{\circ}$', fontsize=18, x=0.1, weight='bold')
    cax = ax.imshow(dt_arr, interpolation='nearest', extent=extent, aspect='auto', vmin=1, vmax=10)

    ax.tick_params(labelsize=16)
    
    if (i <= 11):
        ax.set_xticklabels([])
    else:
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)      
    if (i == 0 or i == 6 or i == 12):
        pass
    else:
        ax.set_yticklabels([])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(cax, cax=cbar_ax)
cb.set_label('Sc$^{\prime}$S delay time (s)', y=0.5, fontsize=26)
cb.ax.tick_params(labelsize=18)

plt.show()
        
fig.savefig('Plots/' + name + '/2d_dt_plots.png', bbox_inches='tight')
fig.savefig('Plots/' + name + '/2d_dt_plots.pdf', bbox_inches='tight')
