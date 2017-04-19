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
    # Initialise counters
    s = 0
    clickCounter = 0

    # Initialise storage for clicks
    delayTimes = []
    amplRatios = []
    azis = []
    dists = []

    # Filter frequencies
    fmin = 0.033
    fmax = 0.1
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.name = StringVar()
        self.event = StringVar()
        self.eventName = StringVar()
        self.azimuth = StringVar()
        self.distance = StringVar()
        ttk.Label(master=root, text="Name").grid(row=1, column=0)
        ttk.Entry(master=root, textvariable=self.name).grid(row=1, column=1, columnspan=2)

        self.createWidgets()

    def initiate(self, canvas, ax):
        self.dir = 'Data/' + str(self.name.get()) + '/'
        self.seislist = glob.glob(self.dir + '/*PICKLE')
        print(self.seislist)

        self.resetCounter()

        self.plot(canvas, ax)

        self.distAzUpdate()

    def counterUpdate(self):
        self.eventName.set(str(self.seislist[self.s]))
        self.event.set(str(self.s) + ' / ' + str(len(self.seislist)))

    def distAzUpdate(self):
        self.azimuth.set('Azimuth : ' + str(self.az))
        self.distance.set('Distance : ' + str(self.dist))
        
    def resetCounter(self):
        print(self.s)
        self.s = 0
        self.clickCounter = 0
        self.counterUpdate()

    def acceptData(self, canvas, ax, event=None):
        self.s += 1
        if (self.s == len(self.seislist)):
            print('Last seismogram')
            self.writeFile()
        else:
            if (self.clickCounter == 0):
                self.delayTimes.append(0)
                self.amplRatios.append(0)
            self.writeFile()
            self.clickCounter = 0
            self.counterUpdate()
            self.plot(canvas, ax)
            self.distAzUpdate()

    def writeFile(self):
        print(self.delayTimes, self.amplRatios)
        peakData = open(self.dir + 'peakData.txt', 'w')
        for i in range(len(self.delayTimes)):
            peakData.write("%s %s %s %s \n" % (self.delayTimes[i], self.amplRatios[i], self.azis[i], self.dists[i]))
        peakData.close()

    def createWidgets(self):
        ttk.Button(master=root, text="Start", command=lambda: self.initiate(canvas,ax), underline=0).grid(row=1, column=3, columnspan=2)        
        ttk.Button(master=root, text="Accept", command=lambda: self.acceptData(canvas,ax)).grid(row=1, column=5, columnspan=2)
        ttk.Button(master=root, text="Reset", command=lambda: self.resetCounter()).grid(row=1, column=7)

        # Create labels so can track which event is being viewed
        ttk.Label(master=root, textvariable=self.event).grid(row=4, column=5)
        ttk.Label(master=root, textvariable=self.eventName).grid(row=4, column=3, columnspan=2)

        # Create labels for distance and azimuth
        ttk.Label(master=root, textvariable=self.azimuth).grid(row=3, column=2, columnspan=2)
        ttk.Label(master=root, textvariable=self.distance).grid(row=3, column=4, columnspan=2)

        # Create canvas
        fig = plt.figure(figsize=(16,8), dpi=100)
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().grid(row=2, column=0, columnspan=8)
        canvas.show()

        self.update()
        
        root.bind("<Right>", lambda _: self.acceptData(canvas, ax))

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

        self.azis.append(self.az)
        self.dists.append(self.dist)

        self.seisT.data = np.gradient(self.seisT.data, self.seisT.stats.delta)

        # Plotting data. Both components normalised by max amplitude on transverse component
        ax.plot(self.seisT.times() + self.tshift, self.seisT.data / np.max(self.seisT.data), 'k', linewidth=2)

        # Plotting travel time predictions
        for k in seis[0].stats.traveltimes.keys():
            if seis[0].stats.traveltimes[k] != None:
                ax.plot(seis[0].stats.traveltimes[k], 0.0, 'g', marker='o', markersize=4)
                ax.text(seis[0].stats.traveltimes[k], -0.2, k, fontsize=8)

        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])

        cid = canvas.mpl_connect('button_press_event', self.__onclick__)
        
        canvas.draw()
        
    def __onclick__(self, click):
        print('Click!')
        self.clickCounter += 1
        self.point = (click.xdata, click.ydata)
        xArr = self.seisT.times() + self.tshift
        if (self.clickCounter == 1):
            self.idx1 = self.findNearest(xArr, click.xdata)
        if (self.clickCounter == 2):
            self.idx2 = self.findNearest(xArr, click.xdata)
            self.calcValues(xArr, self.idx1, self.idx2)            
            
    def findNearest(self, array, value):
        idx = (np.abs(array-value)).argmin() 
        return idx

    # Function that calculates the dt and the amplitude ratio given 2 indices
    def calcValues(self,xArr, idx1, idx2):
        print(xArr[idx2], xArr[idx1])
        delayTime = xArr[idx2] - xArr[idx1]
        amplRatio = self.seisT.data[idx1] / self.seisT.data[idx2]
        print(delayTime, amplRatio)
        self.delayTimes.append(delayTime)
        self.amplRatios.append(amplRatio)
        
root = tk.Tk()
root.title("Part III Project - Parameter Calculator")
app = Application(master=root)
app.mainloop()
