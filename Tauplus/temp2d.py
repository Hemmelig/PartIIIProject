import matplotlib.pyplot as plt
# Obspy is a seismic toolkit
import obspy
print(obspy.__path__)
from obspy.taup import TauPyModel
from obspy.taup.taup import travelTimePlot
from obspy.taup.tau import Arrivals
import obspy.taup
from collections import OrderedDict
# phases.py needs to be in the same directory
import phases


# Phases to plot, e.g. plotphase =["PKJKP", "SKKS"]
plotphase =["S", "ScS","Sdiff", "ScS^2871ScS"]

# depth of earthquake in km
depth_earthquake= 414.89

# Defines array of distances to compute ray path at
# format: np.arange(start, stop, stepsize)
# e.g np.arange(0,4,1)=[0,1,2,3]
dist_raypaths = np.arange(70,100.1, 5)
# Defines array of distances to compute traveltime at
dist_traveltimes = np.arange(70,100.1,1)


# Anomalous layer
top_r = np.arange(2886.0, 2851.0, 1)
bot_r = 2891.0
dv_s = np.arange(-0.1, -0.3, 0.01)
dv_p = -0.1
drho = 0.1

# Color pallette to use for plotting
cm = plt.get_cmap('gist_rainbow')

model_name = 'prem'

# velocity model as a function of depth.
model = obspy.taup.taup_create.TauPCreate(input_filename='prem.nd', output_filename='prem.npz')
model.load_velocity_model()
model.run() # Saves PREM to prem.npz
model.output_filename = 'prem_mod.npz'

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


plt.subplot(1,3,1) #
plt.plot(vs, depth, label='prem') #

model.v_mod.layers = np.array(modified_layers)
model.run()
depth = []
vs = []
for layer in model.v_mod.layers:
    depth.append(layer[0])
    depth.append(layer[1])
    vs.append(layer[4])
    vs.append(layer[5])

plt.plot(vs, depth, label='prem+ULVZ')
plt.gca().invert_yaxis()
plt.legend()

##### plot for PREM in black for reference

model = TauPyModel(model='./prem.npz')

# Loop through phases
for ph, phase in enumerate(plotphase):
    # Loop through distances
    for distance in dist_traveltimes:
        # Print out what is being done
        print('computing time of ' + phase +' at '+ str(distance) + ' deg')
        # Compute ray path and arrival time.
        arrivals = model.get_ray_paths(depth_earthquake, distance, phase_list=[phase])
        # Plot arrival times
        plt.subplot(1,3,2)
        for ind, arr in enumerate(arrivals):
            plt.plot(distance, arr.time,'.',color =cm(ph/float(len(plotphase))), marker ='^', label=phase +'_prem')
        plt.subplot(1,3,3)
        for ind, arr in enumerate(arrivals):
            plt.plot(distance, arr.ray_param,'.',color =cm(ph/float(len(plotphase))), marker ='^',label=phase +'_prem')

##### plot for modified PREM

model = TauPyModel(model='./prem_mod.npz')

# Loop through phases
for ph, phase in enumerate(plotphase):
    # Loop through distances
    for distance in dist_traveltimes:
        # Print out what is being done
        print('computing time of ' + phase +' at '+ str(distance) + ' deg')
        # Compute ray path and arrival time.
        arrivals = model.get_ray_paths(depth_earthquake, distance, phase_list=[phase])
        # Plot arrival times
        plt.subplot(1,3,2)
        for ind, arr in enumerate(arrivals):
                plt.plot(distance, arr.time,'.',color =cm(ph/float(len(plotphase))), marker ='o',label=phase + '_mod')
        plt.subplot(1,3,3)
        for ind, arr in enumerate(arrivals):
            print(arr.time)
            plt.plot(distance, arr.ray_param,'.',color =cm(ph/float(len(plotphase))),marker ='o', label=phase + '_mod')

plt.subplot(1,3,2)
# Annotate travel time figure
plt.xlabel('distance (dg)')
plt.title('distance vs time')
#plt.ylabel('time (s)')
plt.xlim([min(dist_traveltimes), max(dist_traveltimes)])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc=4)
plt.subplot(1,3,3)
# Annotate travel time figure
plt.xlabel('distance (dg)')
#plt.ylabel('time (s)')
plt.ylabel('ray param (s/rad)')
plt.title('distance vs. ray param')

plt.savefig('example.png')
plt.savefig('example.pdf')


# Initialize the figure and axes for the ray path plot
# ax_right is used for paths plotted on the right half.
fig, ax = plt.subplots(1,1,figsize=(8,8), subplot_kw=dict(polar=True))
plt.title('Ray paths in modified model')
ax_right = ax
ax_right.set_theta_zero_location('N')
ax_right.set_theta_direction(-1)
# ax_left is used for paths plotted on the left half.
ax_left = fig.add_axes(ax_right.get_position(), projection='polar',
                       label='twin', frameon=False)
ax_left.set_theta_zero_location('N')
ax_left.set_theta_direction(+1)
ax_left.xaxis.set_visible(False)
ax_left.yaxis.set_visible(False)


# Loop through phases
for ph, phase in enumerate(plotphase):
    # Loop through distances
    for distance in dist_raypaths:
        # Print out what is being done
        print('computing ray path for ' +phase+ ' at '+ str(distance) +' deg')
        # Computer ray path
        arrival = model.get_ray_paths(depth_earthquake, distance, phase_list= [phase])
        # If there is no arrival at this distance, continue to the next, otherwise plot.
        if not len(arrival):
            print('No arrivals at ', distance, ' degrees for ', phase)
            continue
        
        phases.plot(arrival,plot_type='spherical', legend=False, label_arrivals=False,
                    show=False, ax=ax_right, color =cm(ph/float(len(plotphase))))


plt.show()
