import os
import numpy as np
import matplotlib.pyplot as plt
import obspy
print(obspy.__path__)
from obspy.taup import TauPyModel
from obspy.taup.taup import travelTimePlot
from obspy.taup.tau import Arrivals
import obspy.taup
from collections import OrderedDict
import phases

# Phases to plot, e.g. plotphase = ["S", "ScS"]
plotphase = ["S", "ScS"]

# Depth of earthquake in km
depth_earthquake = 414.89

# Defines array of distances to compute ray path at
# format: np.arange(start, stop, stepsize)
# e.g. np.arange(0,4,1)=[0,1,2,3]
dist_raypaths = 90.
print(dist_raypaths)
# Defines array of distances to compute traveltime at
distance = 90.

# Anomalous layer
#top_r_arr = np.arange(2886.0, 2850.0, -1)
top_r_arr = np.arange(2878.0, 2882.0, 1)
bot_r = 2891.0
#dv_s_arr = np.arange(-0.1, -0.26, -0.01)
dv_s_arr = np.arange(-0.24, -0.2, 0.01)
dv_p = -0.1
drho = 0.1

# Colour palette to use for plotting
cm = plt.get_cmap('gist_rainbow')

model_name = 'prem'

print(len(top_r_arr), len(dv_s_arr))
dt_arr = np.zeros((len(top_r_arr), len(dv_s_arr)))

for i in range(len(top_r_arr)):
    top_r = top_r_arr[i]
    
    for j in range(len(dv_s_arr)):
        print(dv_s_arr[j])
        dv_s = dv_s_arr[j]

        
        # Velocity model as a function of depth
        model = obspy.taup.taup_create.TauPCreate(input_filename='prem.nd', output_filename='prem.npz')
        model.load_velocity_model()
        model.run()
        model.output_filename = 'prem_mod_' + str(i) + '_' + str(j) + '.npz'

        depth = []
        vs = []
        modified_layers = []
        set_toplayer = True
        set_botlayer = True

        for layer in model.v_mod.layers:
            depth.append(layer[0])
            depth.append(layer[1])
            vs.append(layer[4])
            vs.append(layer[5])
            if layer[0] >= top_r and layer[1] <= bot_r and set_toplayer == False:
                layer[2] = (1. + dv_p) * layer[2]
                layer[4] = (1. + dv_s) * layer[4]
                layer[6] = (1. + drho) * layer[6]
                layer[3] = (1. + dv_p) * layer[3]
                layer[5] = (1. + dv_s) * layer[5]
                layer[7] = (1. + drho) * layer[7]

            if layer[1] >= top_r and set_toplayer == True:
                layertmp = layer.copy()
                layertmp[1] = top_r
                layer[0] = top_r
                layer[2] = (1. + dv_p) * layer[2]
                layer[4] = (1. + dv_s) * layer[4]
                layer[6] = (1. + drho) * layer[6]
                layer[3] = (1. + dv_p) * layer[3]
                layer[5] = (1. + dv_s) * layer[5]
                layer[7] = (1. + drho) * layer[7]
        
                if layertmp[0] != layertmp[1]:
                    modified_layers.append(layertmp)
                set_toplayer = False

            if layer[1] >= bot_r and set_botlayer == True:
                layertmp = layer.copy()
                layertmp[0] = bot_r
                layer[1] = bot_r
                layertmp[2] = (1. + dv_p) * layer[2]
                layertmp[4] = (1. + dv_s) * layer[4]
                layertmp[6] = (1. + drho) * layer[6]
                layertmp[3] = (1. + dv_p) * layer[3]
                layertmp[5] = (1. + dv_s) * layer[5]
                layertmp[7] = (1. + drho) * layer[7]

                if layertmp[0] != layertmp[1]:
                    modified_layers.append(layertmp)
                set_botlayer = False

            if layer[0] != layer[1]:
                modified_layers.append(layer)

        model.v_mod.layers = np.array(modified_layers)
        model.run()
        depth = []
        vs = []
        for layer in model.v_mod.layers:
            depth.append(layer[0])
            depth.append(layer[1])
            vs.append(layer[4])
            vs.append(layer[5])

        # Plot for modified PREM
        
        model = TauPyModel(model='./prem_mod_' + str(i) + '_' + str(j) + '.npz')
                
        print('Computing travel time of S')
        arrivals_S = model.get_travel_times(depth_earthquake, distance, phase_list=[plotphase[0]])
        for ind, arr in enumerate(arrivals_S):
            S_arr = arr.time
        print(S_arr)
            
        print('Computing travel time of ScS')
        arrivals_ScS = model.get_travel_times(depth_earthquake, distance, phase_list=[plotphase[1]])
        for ind, arr in enumerate(arrivals_ScS):
            ScS_arr = arr.time
            print(ScS_arr)

        dt = ScS_arr - S_arr
        dt_arr[i][j] = dt

        print(i, j)

        del arrivals_S
        del arrivals_ScS

fig = plt.figure(figsize=(6,4))

ax = fig.add_subplot(111)
ax.set_title('colorMap')
#scaler = dt_arr.max()

#def f(x):
#    return x / scaler
#f = np.vectorize(f)

#dt_arr = f(dt_arr)

y, x = np.meshgrid(top_r_arr, dv_s_arr)

print(dt_arr, top_r_arr, dv_s_arr)
plt.pcolormesh(x, y, dt_arr, cmap=cm)
#ax.set_aspect('equal')

plt.colorbar(orientation='vertical')
plt.show()
        
#plt.savefig('example.png')
#plt.savefig('example.pdf')

#plt.show()
