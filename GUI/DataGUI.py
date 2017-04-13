#!/usr/bin/env python3

### The purpose of this code is to provide a simple GUI allowing the user to perform
### a range of different actions on a given dataset
### To use:
### 1. Fill in the details of the event and stations of interest
### 2. Download the data then process it
### 3. Add the CMT source data - open base file, edit then save in relevant Data folder
### 4. Add synthetics and traveltimes
### 5. Plot data. !!! Need to add distance and azimuth functionality of plots !!!

# To do: turn the whole script into a class, better coding standard

#####################################################################################
### Import modules and external python scripts
#####################################################################################

import sys, os, subprocess
from tkinter import *
import tkinter as tk
from tkinter import ttk, filedialog
from download_data_one_event import *
from data_processing import *
from add_synthetics_topickle import *
from add_travel_times import *
from add_event_cmt import *
from plot_data_with_azimuth import *
from plot_data_with_distance import *
from tkinter.scrolledtext import ScrolledText

#####################################################################################
    
root = Tk()
root.title("Part III Project")

# Create menu
menubar = Menu(root)
menubar.add_command(label="Hello!")
menubar.add_command(label="Exit", command=root.quit)
root.config(menu=menubar)

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)


#####################################################################################
### Define all variables being input and create input boxes
#####################################################################################

# Event and station parameters (this is really REALLY ugly..)
name, latitude, longitude, starttime, endtime, maxrad, minmag, maxmag, distmin, distmax, azmin, azmax, lengthoftrace = StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar(), StringVar()

ttk.Label(mainframe, text="Name").grid(column=1, row=1, sticky=W)
name_entry = ttk.Entry(mainframe, width=7, textvariable=name).grid(column=2, row=1, sticky=(W, E), columnspan=2)

# Label names
labels = ["Event Parameters", "Latitude", "Longitude", "Start Time", "End Time", "Max Radius", "Min Magnitude", "Max Magnitude", "Station Parameters", "Min Distance", "Max Distance", "Min Azimuth", "Max Azimuth", "Length of trace"]
parameters = [latitude, longitude, starttime, endtime, maxrad, minmag, maxmag, distmin, distmax, azmin, azmax, lengthoftrace]

rowNo = 3

for label in labels:
    ttk.Label(mainframe, text=label).grid(column=1, row=rowNo, sticky=W)
    rowNo += 1

rowNo = 4

for param in parameters:
    if (rowNo == 6) or (rowNo == 7):
        ttk.Entry(mainframe, width=7, textvariable=param).grid(column=2, row=rowNo, sticky=(W, E), columnspan=2)
    else:
        ttk.Entry(mainframe, width=7, textvariable=param).grid(column=2, row=rowNo, sticky=(W, E))
    if (rowNo == 10):
        rowNo += 2
    else:
        rowNo += 1

del rowNo

def download():
    # Get event parameters
    n = str(name.get())
    l1 = float(latitude.get())
    l2 = float(longitude.get())
    st = str(starttime.get())
    et = str(endtime.get())
    mr = float(maxrad.get())
    mim = float(minmag.get())
    mam = float(maxmag.get())

    # Get station parameters
    dmin = float(distmin.get())
    dmax = float(distmax.get())
    amin = float(azmin.get())
    amax = float(azmax.get())
    ltrace = float(lengthoftrace.get())

    # Pass arguments to download data script
    downloadData(n, l1, l2, st, et, mr, mim, mam, dmin, dmax, amin, amax, ltrace)

def process():
    # Get event name
    n = str(name.get())

    # Pass name to data processing script
    processData(n)

# Buttons to execute download of data and data processing
ttk.Button(mainframe, text="Download", command=download).grid(column=3, row=15, sticky=W)
ttk.Button(mainframe, text="Process", command=process).grid(column=3, row=16, sticky=W)   

#####################################################################################


#####################################################################################
### Section for adding synthetics
#####################################################################################

# Create frame
synthframe = ttk.Frame(mainframe, padding="12 12 12 12", relief=RAISED)
synthframe.grid(column=1, row=17, columnspan=2, rowspan=3, sticky=(N, W, E, S))
synthframe.columnconfigure(0, weight=1)
synthframe.rowconfigure(0, weight=1)

def addSynths():
    # Get event name
    n = str(name.get())
    clean = bool(v.get())
    
    addSynthetics(n, clean)
    
# Create label, radiobuttons and execute button
v = BooleanVar()
ttk.Label(synthframe, text="Remove synthetics?").pack(anchor=W)
Radiobutton(synthframe, text="Yes", variable=v, value=True).pack(anchor=W)
Radiobutton(synthframe, text="No", variable=v, value=False).pack(anchor=W)
ttk.Button(synthframe, text="Add Synthetics", command=addSynths).pack(anchor=W)

#####################################################################################


#####################################################################################
### Section for adding travel times
#####################################################################################

# Create frame
travelframe = ttk.Frame(mainframe, padding="12 12 12 12", relief=RAISED)
travelframe.grid(column=4, row=12, columnspan=2, rowspan=5, sticky=(N, W, E, S))
travelframe.columnconfigure(0, weight=1)
travelframe.rowconfigure(0, weight=1)

def read(phase):
    return vPhases[phase].get()

def addTravelTimes():
    # Get event name
    n = str(name.get())

    # Check for any phases to include
    phasesToInclude = []
    for i in range(len(vPhases)):
        if read(i):
            phasesToInclude.append(phases[i])
        else:
            continue

    travelTimes(n, phasesToInclude)
    
# Create label, checkbuttons and execute button
ttk.Label(travelframe, text="Travel times").pack(anchor=W)
ttk.Label(travelframe, text="Select phases").pack(anchor=W)

# Loop through phase options and create check buttons for each. Just add more names
# to the phases list to produce more buttons. May need to add a counter and if
# statement in order to split into two columns beyond a certain number of phases
phases = ["S", "Sdiff", "ScS", "SKS", "SKKS"]
vPhases = [None] * len(phases)

for i in range(len(phases)):
    vPhases[i] = BooleanVar()
    Checkbutton(travelframe, text=phases[i], variable=vPhases[i]).pack(anchor=W)

ttk.Button(travelframe, text="Add Travel Times", command=addTravelTimes).pack(anchor=W)

#####################################################################################


#####################################################################################
### Section to add the CMT to the data
#####################################################################################

# Create frame
cmtframe = ttk.Frame(mainframe, padding="12 12 12 12")
cmtframe.grid(column=4, row=1, rowspan=11)
cmtframe.columnconfigure(0, weight=1)
cmtframe.rowconfigure(0, weight=1)

def open_command():
    myfile = filedialog.askopenfile(parent=mainframe,
                              title='Open a Python file', mode='r')
    loadedfile = myfile.read()
    myfile.close()
    cmt_entry.insert("end", loadedfile)

def save_command():
    file = filedialog.asksaveasfile(mode='w')
    if file != None:
    # Slice off the last character from get, as an extra return is added
        data = cmt_entry.get('1.0', END+'-1c')
        file.write(data)
        file.close()

def addCMT():
    # Get event name
    n = str(name.get())   

    addEventCMT(n)
    
ttk.Label(cmtframe, text="CMT Source Data").pack(anchor=N)

cmt_entry = ScrolledText(cmtframe)
cmt_entry.pack(anchor=N)

cmtsaveframe = ttk.Frame(cmtframe, padding="12 12 12 12")
cmtsaveframe.pack()

ttk.Button(cmtsaveframe, text="Open", command=open_command).pack(anchor=W)
ttk.Button(cmtsaveframe, text="Save", command=save_command).pack(anchor=W)
ttk.Button(cmtsaveframe, text="Add", command=addCMT).pack(anchor=E)

#####################################################################################

#####################################################################################
### Section to add the CMT to the data
#####################################################################################

### Want to extend this section to give better choices: amplitude, distance divisions
### for azimuth plotter etc

def plotter():
    # Get event name
    n = str(name.get())
    if plotType.get():
        plotWithAzimuth(n)
    else:
        plotWithDistance(n)

plotType = BooleanVar()
Radiobutton(mainframe, text="Azimuth", variable=plotType, value=True).grid(column=4, row=17, columnspan=2, sticky=W)
Radiobutton(mainframe, text="Distance", variable=plotType, value=False).grid(column=4, row=18, columnspan=2, sticky=W)

ttk.Button(mainframe, text="Plot", command=plotter).grid(column=5, row=18, columnspan=2, sticky=W)

#####################################################################################

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

#name_entry.focus()

root.mainloop()
