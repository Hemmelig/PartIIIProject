import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import csv

name = sys.argv[1]

# Defines array of distances to compute traveltime at
distance = sys.argv[2]

# Anomalous layer
top_r_arr = np.arange(2856.0, 2886.5, 0.5)
dv_s_arr = np.arange(-0.24, -0.14, 0.002)

dir = 'Data/PeakData/' + name

with open(dir + '_' + str(int(distance)) + '_2d_dt.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    dt_arr = [[float(e) for e in r] for r in reader]

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.set_title('Delay time against velocity reduction and height variation')
cax = ax.matshow(dt_arr)
fig.colorbar(cax)

# Convert array values into strings
top_r_arr = [str(x) for x in top_r_arr]
dv_s_arr = [str(x) for x in dv_s_arr]

ax.set_xticklabels([''] + dv_s_arr)
ax.set_yticklabels([''] + top_r_arr)

plt.show()
        
fig.savefig('Plots/' + name + '/2d_dt_plot_' + str(int(distance)) + '.png')
fig.savefig('Plots/' + name + '/2d_dt_plot_' + str(int(distance)) + '.pdf')
