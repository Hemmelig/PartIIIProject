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
from tkinter.scrolledtext import ScrolledText


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
    clickCounter = 0

    # Initialise storage for clicks
    pSs = []
    pScSs = []

    # Filter frequencies
    fmin = 0.
    fmax = 0.
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.name = StringVar()
        self.event = StringVar()
        self.eventName = StringVar()
        self.pScSav = StringVar()
        self.pSav = StringVar()
        self.pCount = StringVar()
        self.stdScS = StringVar()
        self.stdS = StringVar()
        self.distance = StringVar()
        self.azimuth = StringVar()
        self.freqBand = IntVar()
        self.freqBand.set(1)
        self.freqMin = DoubleVar()
        self.freqMax = DoubleVar()
        self.dataOrSyn = IntVar()
        self.dataOrSyn.set(1)

        self.createWidgets()

    def createWidgets(self):
        # Create menu
        menubar = Menu(root)
        menubar.add_command(label="Instructions", command=lambda: self.instructions())
        menubar.add_command(label="Exit", command=root.quit)
        root.config(menu=menubar)
        
        # Create options section
        optionFrame = ttk.Frame(root, padding="12 12 12 12", relief=RAISED)
        optionFrame.grid(column=0, row=2, columnspan=2, rowspan=10, sticky=(N, W, E, S))
        optionFrame.columnconfigure(0, weight=1)
        optionFrame.rowconfigure(0, weight=1)

        # Option label
        ttk.Label(master=optionFrame, text="Options").pack(anchor=W)
        
        # Name label and entry box
        ttk.Label(master=optionFrame, text="Name").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.name).pack(anchor=W)
        
        # Real data or synthetics option toggle
        Radiobutton(master=optionFrame, text="Real data", variable=self.dataOrSyn, value=1).pack(anchor=W)
        Radiobutton(master=optionFrame, text="Synthetics", variable=self.dataOrSyn, value=2).pack(anchor=W)

        # Control buttons for plots
        ttk.Button(master=optionFrame, text="Start", command=lambda: self.initiate(canvas,ax)).pack(anchor=W)       
        ttk.Button(master=optionFrame, text="Accept", command=lambda: self.acceptData(canvas,ax)).pack(anchor=W)
        ttk.Button(master=optionFrame, text="Reset", command=lambda: self.resetCounter()).pack(anchor=W)

        # Frequency options and update button
        Radiobutton(master=optionFrame, text="10-30s", variable=self.freqBand, value=1).pack(anchor=W)
        Radiobutton(master=optionFrame, text="10-20s", variable=self.freqBand, value=2).pack(anchor=W)
        Radiobutton(master=optionFrame, text="Custom", variable=self.freqBand, value=3).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Min. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMax).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Max. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMin).pack(anchor=W)

        ttk.Button(master=optionFrame, text="Update", command=lambda: self.updateFreqBand(canvas, ax)).pack(anchor=W)

        # Create labels for distance and azimuth
        ttk.Label(master=optionFrame, text="Azimuth").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.azimuth).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Distance").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.distance).pack(anchor=W)

        # Create label for current period average
        ttk.Label(master=optionFrame, text="Samples").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.pCount).pack(anchor=W)
        
        # Create label for current period average
        ttk.Label(master=optionFrame, text="ScS period").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.pScSav).pack(anchor=W)

        # Create label for current period average
        ttk.Label(master=optionFrame, text="ScS std").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.stdScS).pack(anchor=W)
        
        # Create label for current period average
        ttk.Label(master=optionFrame, text="S period").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.pSav).pack(anchor=W)

        # Create label for current period average
        ttk.Label(master=optionFrame, text="S std").pack(anchor=W)
        ttk.Label(master=optionFrame, textvariable=self.stdS).pack(anchor=W) 

        # Add some padding to all widgets in this frame
        for child in optionFrame.winfo_children():
            child.pack_configure(padx=5, pady=5)

        # Create real data canvas
        fig = plt.figure(figsize=(16,8), dpi=100)
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().grid(row=2, column=2, rowspan=6, columnspan=8)
        canvas.show()
        
        self.update()
        
        root.bind("<Right>", lambda _: self.acceptData(canvas, ax))        

    def initiate(self, canvas, ax):
        self.dir = 'Data/' + str(self.name.get()) + '/'
        self.seislist = glob.glob(self.dir + '/*PICKLE')
        print(self.seislist)

        self.resetCounter()
        
        self.updateFreqBand(canvas, ax)

        self.distAzUpdate()

    def updateFreqBand(self, canvas, ax):
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

        if (self.dataOrSyn.get() == 1):
            self.plot(canvas, ax)
        if (self.dataOrSyn.get() == 2):
            self.plotCSEM(canvas, ax)

    def counterUpdate(self):
        self.eventName.set(str(self.seislist[self.s]))
        self.event.set(str(self.s + 1) + ' / ' + str(len(self.seislist)))

    def distAzUpdate(self):
        self.azimuth.set(str(round(self.az, 3)))
        self.distance.set(str(round(self.dist, 3)))
        self.pScSav.set(str(2 * np.mean(self.pScSs)))
        self.pSav.set(str(2 * np.mean(self.pSs)))
        self.pCount.set(str(len(self.pScSs)))
        self.stdS.set(str(2 * np.std(self.pSs)))
        self.stdScS.set(str(2 * np.std(self.pScSs)))
        
    def resetCounter(self):
        self.s = 0
        self.clickCounter = 0
        self.counterUpdate()

    def acceptData(self, canvas, ax, event=None):
        self.s += 1
        
        if ((self.s) == len(self.seislist)):
            print('Mean S period: ' + str((2 * np.mean(self.pSs))))
            print('Standard deviation: ' + str(np.std(self.pSs)))
            print('Mean ScS period: ' +  str((2 * np.mean(self.pScSs))))
            print('Standard deviation: ' + str(np.std(self.pScSs)))
            print('Last seismogram')
        else:
            print(self.s + 1, len(self.seislist))
            self.clickCounter = 0
            self.counterUpdate()
            self.updateFreqBand(canvas, ax)
            self.distAzUpdate()

    def plot(self, canvas, ax):
        # Clear current plot
        ax.clear()

        # Get current data
        seis = read(self.seislist[int(self.s)], format='PICKLE')
        
        seis.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
        seis.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)
        
        phase = 'Sdiff'
        if seis[0].stats.traveltimes[phase] is None:
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

        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])

        ax.set_title('%s' % (str(self.eventName.get())), loc='left')
        ax.set_title('%s' % (str(self.event.get())), loc='right')

        cid = canvas.mpl_connect('button_press_event', self.__onclick__)
        
        canvas.draw()

    def plotCSEM(self, canvas, ax):
        # Clear current plot
        ax.clear()

        # Get current data
        seis = read(self.seislist[int(self.s)], format='PICKLE')
        sta = seis[0].stats.station
        nw = seis[0].stats.network
        slat = seis[0].stats['stla']
        slon = seis[0].stats['stlo']
        self.tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        self.az = seis[0].stats['az']
        self.dist = seis[0].stats['dist']

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
            self.times.append(float(columns[0]) - seis[0].stats.traveltimes[phase])
            self.seistoplot.append(float(columns[1]))     

        norm = np.max(np.absolute(self.seistoplot))
        self.seistoplot = [x / norm for x in self.seistoplot]
            
        ax.plot(self.times, self.seistoplot, 'k', linewidth=2) 

        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([430., 510.])

        ax.set_title('%s' % (str(self.eventName.get())), loc='left')
        ax.set_title('%s' % (str(self.event.get())), loc='right')

        cid = canvas.mpl_connect('button_press_event', self.__onclick__)

        canvas.draw()
        
    def __onclick__(self, click):
        print('Click!')
        self.clickCounter += 1
        self.point = (click.xdata, click.ydata)
        if (self.dataOrSyn.get() == 1):
            xArr = self.seisT.times() + self.tshift
            array = self.seisT.data
        if (self.dataOrSyn.get() == 2):
            xArr = self.times
            array = self.seistoplot
        print(click.xdata)
            
        if (self.clickCounter == 1):
            idxtmp = (np.abs(xArr - click.xdata)).argmin()
            self.idx1 = self.findMin(array, idxtmp)
        elif (self.clickCounter == 2):
            idxtmp = (np.abs(xArr - click.xdata)).argmin()
            self.idx2 = self.findMax(array, idxtmp)
        elif (self.clickCounter == 3):
            idxtmp = (np.abs(xArr - click.xdata)).argmin()
            self.idx3 = self.findMin(array, idxtmp)
        elif (self.clickCounter == 4):
            idxtmp = (np.abs(xArr - click.xdata)).argmin()
            self.idx4 = self.findMax(array, idxtmp)
            self.calcValues(xArr, self.idx1, self.idx2, self.idx3, self.idx4)
            
    def findMax(self, array, idxtmp):
        while (np.abs(array[idxtmp]) < np.abs(array[idxtmp + 1])):
            idxtmp += 1
        while (np.abs(array[idxtmp]) < np.abs(array[idxtmp - 1])):
            idxtmp -= 1
        return idxtmp

    def findMin(self, array, idxtmp):
        while (array[idxtmp] > array[idxtmp + 1]):
            idxtmp += 1
        while (array[idxtmp] > array[idxtmp - 1]):
            idxtmp -= 1
        return idxtmp                      

    # Function that calculates the dt and the amplitude ratio given 2 indices
    def calcValues(self, xArr, idx1, idx2, idx3, idx4):
        pS = xArr[idx2] - xArr[idx1]
        pScS = xArr[idx4] - xArr[idx3]

        self.pSs.append(pS)
        self.pScSs.append(pScS)
        
    def instructions(self):
        toplevel = Toplevel()
        label1 = Label(toplevel, text=self.ABOUT_TEXT, height=0, width=100)
        label1.pack(anchor=W)
        label2 = Label(toplevel, text=self.DISCLAIMER, height=0, width=100)
        label2.pack()

root = tk.Tk()
root.title("Part III Project - Period Calculator")
app = Application(master=root)
app.mainloop()
