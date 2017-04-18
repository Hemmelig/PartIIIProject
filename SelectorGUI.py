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
    # Initialise counter
    s = 0

    # Plot synthetics?
    syn = True

    # Filter frequencies
    fmin = 0.033
    fmax = 0.1

    replot = False
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.name = StringVar()
        self.event = StringVar()
        self.eventName = StringVar()
        ttk.Label(master=root, text="Name").grid(row=1, column=0)
        nameEntry = ttk.Entry(master=root, textvariable=self.name)
        nameEntry.grid(row=1, column=1, columnspan=2)
        nameEntry.focus_set()

        self.createWidgets()

    def initiate(self, canvas, ax, canvas2, ax2):
        dir = 'Data/' + str(self.name.get()) + '/'
        self.seislist = glob.glob(dir + '/*PICKLE')
        print(self.seislist)

        self.resetCounter()
        
        # Create directory to dump data into
        self.dirdump = dir + 'Dump'
        if not os.path.exists(self.dirdump):
            os.makedirs(self.dirdump)

        self.plot(canvas, ax, canvas2, ax2)

    def counterUpdate(self):
        self.eventName.set(str(self.seislist[self.s]))
        self.event.set(str(self.s) + ' / ' + str(len(self.seislist)))
        
    def resetCounter(self):
        print(self.s)
        self.s = 0
        self.clickCounter = 0
        self.counterUpdate()

    def rejectData(self, canvas, ax, canvas2, ax2, event=None):
        shutil.move(self.seislist[self.s], self.dirdump)
        self.s += 1
        self.clickCounter = 0
        self.counterUpdate
        self.plot(canvas, ax, canvas2, ax2)

    def acceptData(self, canvas, ax, canvas2, ax2, event=None):
        self.s += 1
        self.clickCounter = 0
        self.counterUpdate
        print('Data accepted')
        self.plot(canvas, ax, canvas2, ax2)

    def invertData(self, canvas, ax, canvas2, ax2, event=None):
        self.seisR.data = self.seisR.data * -1.
        self.seisT.data = self.seisT.data * -1.
              
        if len(self.seisZ) > 1:
            self.seisZ = self.seisZ.select(location='10')
            self.seisZ[0].stats = self.seisStats
            
        seisnew = Stream()
        seisnew.append(self.seisZ[0])
        seisnew.append(self.seisR)
        seisnew.append(self.seisT)
        seisnew.append(self.synT)
        seisnew.append(self.synR)
        seisnew.append(self.synZ)
        # Write out in PICKLE format
        filename = self.seislist[self.s]
        seisnew.write(filename, 'PICKLE')
        self.replot = True
        self.plot(canvas, ax, canvas2, ax2)

    def createWidgets(self):
        ttk.Button(master=root, text="Start", command=lambda: self.initiate(canvas,ax,canvas2,ax2)).grid(row=1, column=3)        
        ttk.Button(master=root, text="Accept", command=lambda: self.acceptData(canvas,ax,canvas2,ax2)).grid(row=1, column=4)
        ttk.Button(master=root, text="Reject", command=lambda: self.rejectData(canvas,ax,canvas2,ax2)).grid(row=1, column=5)
        ttk.Button(master=root, text="Invert", command=lambda: self.invertData(canvas,ax,canvas2,ax2)).grid(row=1, column=6)
        ttk.Button(master=root, text="Reset", command=lambda: self.resetCounter()).grid(row=1, column=7)

        # Create labels so can track which event is being viewed
        ttk.Label(master=root, textvariable=self.event).grid(row=4, column=5)
        ttk.Label(master=root, textvariable=self.eventName).grid(row=4, column=3, columnspan=2)

        # Create canvas
        fig = plt.figure(figsize=(8,4), dpi=100)
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().grid(row=2, column=0, columnspan=8)
        canvas.show()

        fig2 = plt.figure(figsize=(8,4), dpi=100)
        ax2 = fig2.add_axes([0.1,0.1,0.8,0.8])
        canvas2 = FigureCanvasTkAgg(fig2, master=root)
        canvas2.get_tk_widget().grid(row=3, column=0, columnspan=8)
        canvas2.show()

        self.update()
        
        root.bind("<Right>", lambda _: self.acceptData(canvas, ax, canvas2, ax2))
        root.bind("<Left>", lambda _: self.rejectData(canvas, ax, canvas2, ax2))
        root.bind("<Up>", lambda _: self.invertData(canvas, ax, canvas2, ax2))      

    def plot(self, canvas, ax, canvas2, ax2):
        # Clear current plot
        ax.clear()
        ax2.clear()

        # Get current data
        print(self.s, len(self.seislist))
        seis = read(self.seislist[int(self.s)], format='PICKLE')

        if (self.replot == False):
            seis.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
            seis.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)
        phase = 'Sdiff'
        if seis[0].stats.traveltimes[phase] is None:
            phase = 'S'

        tshift = seis[2].stats['starttime'] - seis[2].stats['eventtime']

        self.seisR = seis.select(channel='BHR')[0]
        self.seisT = seis.select(channel='BHT')[0]
        self.seisZ = seis.select(channel='BHZ')
        self.seisStats = seis[0].stats

        # Differentiate data
        self.seisT.data = np.gradient(self.seisT.data, self.seisT.stats.delta)
        self.seisR.data = np.gradient(self.seisR.data, self.seisT.stats.delta)
        

        # Plotting data. Both components normalised by max amplitude on transverse component
        ax.plot(self.seisT.times() + tshift, self.seisT.data / np.max(self.seisT.data), 'k', linewidth=2)
        ax2.plot(self.seisR.times() + tshift, self.seisR.data / np.max(self.seisT.data), 'k', linewidth=2)

        # Retrieve synthetic data
        self.synT = seis.select(channel='BXT')[0]
        self.synR = seis.select(channel='BXR')[0]
        self.synZ = seis.select(channel='BXZ')[0]
        
        # Plot synthetic data
        if (self.syn == True):
            if (self.replot == False):
                self.synT.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
                self.synT.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)

                self.synR.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
                self.synR.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)

            self.synR.data = np.gradient(self.synR.data, self.synR.stats.delta)
            self.synT.data = np.gradient(self.synT.data, self.synT.stats.delta)
                
            ax2.plot(self.synR.times(), self.synR.data / np.max(self.seisT.data), color=[0.5,0.5,0.5])
            ax.plot(self.synT.times(), self.synT.data / np.max(self.seisT.data), color=[0.5,0.5,0.5])

        # Plotting travel time predictions
        for k in seis[0].stats.traveltimes.keys():
            if seis[0].stats.traveltimes[k] != None:
                ax.plot(seis[0].stats.traveltimes[k], 0.0, 'g', marker='o', markersize=4)
                ax.text(seis[0].stats.traveltimes[k], -0.2, k, fontsize=8)

        #ax.title(str(self.seis[0].stats['dist']) + '   ' + str(self.seis[0].stats['az']))
        ax.set_ylim([-1.0, 1.0])
        ax.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])

        ax2.set_ylim([-1.0, 1.0])
        ax2.set_xlim([seis[0].stats.traveltimes[phase] - 50, seis[0].stats.traveltimes[phase] + 100])
        
        canvas.draw()
        canvas2.draw()

        self.replot = False           
    
root = tk.Tk()
root.title("Part III Project - Data Selector")
app = Application(master=root)
app.mainloop()
