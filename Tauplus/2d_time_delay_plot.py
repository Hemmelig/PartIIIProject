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
import csv
import emailer
import sys

# Set email to send notification to.
email = "conor.bacon@gmail.com"

name = sys.argv[1]

# Phases to plot, e.g. plotphase = ["S", "ScS"]
plotphase = ["S", "ScS"]

# Depth of earthquake in km
depth_earthquake = 414.89

# Defines array of distances to compute ray path at
# format: np.arange(start, stop, stepsize)
# e.g. np.arange(0,4,1)=[0,1,2,3]
dist_raypaths = 87.
print(dist_raypaths)
# Defines array of distances to compute traveltime at
distance = 87.

height = 35
maxdv = -0.24
mindv = -0.14
heightStep = 0.5
dvStep = 0.002

# Anomalous layer
top_r_arr = np.arange((2890.5 - height), 2886.0, heightStep)
bot_r = 2891.0
dv_s_arr = np.arange(maxdv, mindv, dvStep)
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
        model = obspy.taup.taup_create.TauPCreate(input_filename='prem.nd', output_filename='./prem/prem.npz')
        model.load_velocity_model()
        model.run()
        model.output_filename = './prem/' + name + '_' + str(int(distance)) + '_prem_mod_' + str(i) + '_' + str(j) + '.npz'

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
        
        model = TauPyModel(model='./prem/' + name + '_' + str(int(distance)) + '_prem_mod_' + str(i) + '_' + str(j) + '.npz')
                
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

dir = '../Data/PeakData/'

with open(dir + name + '_' + str(int(distance)) + '_2d_dt.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in dt_arr]

message = """Your code has finished computing the traveltime grid for:
    Depth Range: 5 - %s km
    (Steps: %s km)
    
    Velocity reduction: %s - %s %%
    (Steps: %s %%)
    
    Distance: %s degrees
    """ % (str(height), str(heightStep), str(maxdv * 100), str(mindv * 100), str(dvStep * 100), str(distance))
    
emailer.sendEmail(email, message)
    
plt.show()
