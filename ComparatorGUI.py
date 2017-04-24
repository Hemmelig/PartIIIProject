#!/usr/bin/env python3

import obspy
from obspy import read, UTCDateTime
from obspy.core import Stream
import os
import os.path
import time
import glob
import shutil
import sys
import subprocess
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import *


'''
This script can be used to select data.
It will loop through the seismograms and show the Radial component (red) and Transverse component (black).
Press the buttons in  the GUI to process the data:
1. Keep data (also bound to left arrow key)
2. Discard data (also bound to right arrow key)
3. Invert polarity (also bound to the up arrow key)
Discarded data will be moved to a directory called Dump. It might be wise to implement a method of retrieval that can be incorporated into the GUI.
'''

class Application(tk.Frame):
    ABOUT_TEXT = """Instructions for use:

    This software can be used to download and process data from the IRIS
    catalogue website. """
    
    DISCLAIMER = """Author

    This interface was developed by Me."""
    
    # Initialise counters
    s = 0

    # Filter frequencies
    fmin = 0.
    fmax = 0.
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.name = StringVar()
        self.event = StringVar()
        self.eventName = StringVar()
        self.azimuth = StringVar()
        self.distance = StringVar()
        self.freqBand = IntVar()
        self.freqBand.set(1)
        self.freqMin = DoubleVar()
        self.freqMax = DoubleVar()

        self.createWidgets()

    def createWidgets(self):
        # Create menu
        menubar = Menu(root)
        menubar.add_command(label="Instructions", command=lambda: self.instructions())
        menubar.add_command(label="Exit", command=root.quit)
        root.config(menu=menubar)
        
        # Create options section
        optionFrame = ttk.Frame(root, padding="12 12 12 12", relief=RAISED)
        optionFrame.grid(column=0, row=2, columnspan=2, rowspan=12, sticky=(N, W, E, S))
        optionFrame.columnconfigure(0, weight=1)
        optionFrame.rowconfigure(0, weight=1)

        # Option label
        ttk.Label(master=optionFrame, text="Options").pack(anchor=W)
        
        # Name label and entry box
        ttk.Label(master=optionFrame, text="Name").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.name).pack(anchor=W)

        # Control buttons for plots
        ttk.Button(master=optionFrame, text="Start", command=lambda: self.initiate(canvas,ax,canvas2,ax2)).pack(anchor=W)       
        ttk.Button(master=optionFrame, text="Accept", command=lambda: self.acceptData(canvas,ax,canvas2,ax2)).pack(anchor=W)
        ttk.Button(master=optionFrame, text="Reset", command=lambda: self.resetCounter()).pack(anchor=W)

        # Frequency options and update button
        Radiobutton(master=optionFrame, text="10-30s", variable=self.freqBand, value=1).pack(anchor=W)
        Radiobutton(master=optionFrame, text="10-20s", variable=self.freqBand, value=2).pack(anchor=W)
        Radiobutton(master=optionFrame, text="Custom", variable=self.freqBand, value=3).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Min. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMax).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Max. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMin).pack(anchor=W)

        ttk.Button(master=optionFrame, text="Update", command=lambda: self.updateFreqBand(canvas, ax, canvas2, ax2)).pack(anchor=W)

        # Create labels for distance and azimuth
        ttk.Label(master=optionFrame, text="Azimuth").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.azimuth).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Distance").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.distance).pack(anchor=W)

        # Create button to make special note for a particular seismogram
        ttk.Button(master=optionFrame, text="Interest", command=lambda: self.addInterest()).pack(anchor=W)
        
        # Add some padding to all widgets in this frame
        for child in optionFrame.winfo_children():
            child.pack_configure(padx=5, pady=5)

        # Create real data canvas
        fig = plt.figure(figsize=(15,5), dpi=100)
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().grid(row=2, column=2, rowspan=6, columnspan=8)
        canvas.show()
        
        # Create synthetics canvas
        fig2 = plt.figure(figsize=(15,5), dpi=100)
        ax2 = fig.add_axes([0.1,0.1,0.8,0.8])
        canvas2 = FigureCanvasTkAgg(fig, master=root)
        canvas2.get_tk_widget().grid(row=8, column=2, rowspan=6, columnspan=8)
        canvas2.show()
        
        self.update()
        
        root.bind("<Right>", lambda _: self.acceptData(canvas, ax, canvas2, ax2))

    def initiate(self, canvas, ax, canvas2, ax2):
        self.dir = 'Data/' + str(self.name.get()) + '/'
        self.seislist = glob.glob(self.dir + '/*PICKLE')
        print(self.seislist)

        self.resetCounter()
        
        self.updateFreqBand(canvas, ax, canvas2, ax2)

        self.distAzUpdate()

    def updateFreqBand(self, canvas, ax, canvas2, ax2):
        freqBand = self.freqBand.get()
        if (freqBand == 1):
            self.fmin = 0.033
            self.fmax = 0.1
        if (freqBand == 2):
            self.fmin = 0.05
            self.fmax = 0.1
        if (freqBand == 3):
            self.fmin = 1 / self.freqMin.get()
            self.fmax = 1 / self.freqMax.get()

        self.plot(canvas, ax)
        self.plotCSEM(canvas2, ax2)

    def counterUpdate(self):
        self.eventName.set(str(self.seislist[self.s]))
        self.event.set(str(self.s + 1) + ' / ' + str(len(self.seislist)))

    def distAzUpdate(self):
        self.azimuth.set(str(round(self.az, 3)))
        self.distance.set(str(round(self.dist, 3)))
        
    def resetCounter(self):
        self.s = 0
        self.counterUpdate()

    def acceptData(self, canvas, ax, canvas2, ax2, event=None):
        self.s += 1
        
        if ((self.s) == len(self.seislist)):
            print('Last seismogram')
        else:
            self.counterUpdate()
            self.updateFreqBand(canvas, ax, canvas2, ax2)
            self.distAzUpdate()
                
    def plot(self, canvas, ax):
        # Clear current plot
        ax.clear()

        # Get current data
        seis = read(self.seislist[int(self.s)], format='PICKLE')
        
        seis.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
        seis.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)
        
        phase = 'S'

        self.tshift = seis[2].stats['starttime'] - seis[2].stats['eventtime']

        self.seisT = seis.select(channel='BHT')[0]

        self.az = seis[0].stats['az']
        self.dist = seis[0].stats['dist']

        self.seisT.data = np.gradient(self.seisT.data, self.seisT.stats.delta)

        # Plotting data. Both components normalised by max amplitude on transverse component
        ax.plot(self.seisT.times() + self.tshift, self.seisT.data / np.max(abs(self.seisT.data)), 'k', linewidth=2)

        # Plotting travel time predictions
        for k in seis[0].stats.traveltimes.keys():
            if seis[0].stats.traveltimes[k] != None:
                ax.plot(seis[0].stats.traveltimes[k], 0.0, 'g', marker='o', markersize=4)
                ax.text(seis[0].stats.traveltimes[k], -0.2, k, fontsize=8)
                
        ax.set_title('%s' % (str(self.eventName.get())), loc='left')
        ax.set_title('%s' % (str(self.event.get())), loc='right')
        
        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])
        
        canvas.draw()

    def plotCSEM(self, canvas2, ax2):
        # Clear current plot
        ax2.clear()

        # Get current data
        seis = read(self.seislist[int(self.s)], format='PICKLE')
        sta = seis[0].stats.station
        nw = seis[0].stats.network
        slat = seis[0].stats['stla']
        slon = seis[0].stats['stlo']
        self.tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        while len(nw) < 2:
            nw = '_' + nw
        while len(sta) < 4:
            sta = '_' + sta
        if len(sta) > 4:
            sta = sta[:4]

        station_name = self.dir + 'Synthetics/UT_%s_%s' % (nw, sta)

        phase = "S"
        self.times = []
        self.seistoplot = []
        
        crs = open(station_name, "r")
        for columns in (row.strip().split() for row in crs):
            self.times.append(float(columns[0]) - 450)
            self.seistoplot.append(float(columns[1]))     

        norm = np.max(np.absolute(self.seistoplot))
        self.seistoplot = [x / norm for x in self.seistoplot]
        #self.times = [x + self.tshift for x in self.times]
        #self.times = [x + seis[0].stats.traveltimes[phase] for x in self.times]
            
        ax2.plot(self.times, self.seistoplot, 'k', linewidth=2) 

        ax2.set_ylim([-1.0, 1.0])
        ax2.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])        

        canvas2.draw()

    def instructions(self):
        toplevel = Toplevel()
        label1 = Label(toplevel, text=self.ABOUT_TEXT, height=0, width=100)
        label1.pack(anchor=W)
        label2 = Label(toplevel, text=self.DISCLAIMER, height=0, width=100)
        label2.pack()        
        
root = tk.Tk()
root.title("Part III Project - Data and Synthetics Comparator")
app = Application(master=root)
app.mainloop()
