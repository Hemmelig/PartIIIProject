#!/usr/bin/env python3

#-----------------------------------------------------------------------------------#
### Import modules and external python scripts -------------------------------------#
#-----------------------------------------------------------------------------------#

import sys, os, subprocess
import obspy
from obspy import read, UTCDateTime
from obspy.core import Stream
import os.path
import time
import glob
import shutil
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
import tkinter as tk
from tkinter import ttk, filedialog
from tkinter.scrolledtext import ScrolledText
sys.path.insert(0, '/raid2/cb793/GUI')
from deconvolve_source import *


#-----------------------------------------------------------------------------------#

class Application(tk.Frame):
    ABOUT_TEXT = """Instructions for use:

    This software can be used to select a source function
    and deconvolve datasets."""
    
    DISCLAIMER = "Author: Conor Bacon"
    
    # Initialise counter
    clickCounter = 0

    # Initialise storage for clicks
    sourceWindow = False

    # Filter frequencies
    fmin = 0.
    fmax = 0.    

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        # Initialise variables
        self.name = StringVar()
        self.dataType = IntVar()
        self.freqBand = IntVar()
        self.deconType = IntVar()

        self.dataType.set(1)
        self.freqBand.set(1)
        self.deconType.set(1)

        self.value = IntVar()
        self.clickCounter = IntVar()

        self.freqMin = DoubleVar()
        self.freqMax = DoubleVar()

        self.waterLevel = DoubleVar()
        self.maxBumps = DoubleVar()
        
        self.createWidgets()

    def createWidgets(self):
        # Create menu
        menubar = Menu(root)
        menubar.add_command(label="Instructions", command=lambda: self.instructions())
        menubar.add_command(label="Exit", command=root.quit)
        root.config(menu=menubar)

        # Create frames and populate with widgets
        # Mainframe
        mainFrame = ttk.Frame(root, padding="3 3 3 3")
        mainFrame.grid(column=0, row=0, sticky=(N, W, E, S))
        mainFrame.columnconfigure(0, weight=1)
        mainFrame.rowconfigure(0, weight=1)

        # Name and processing button frame
        nameFrame = ttk.Frame(mainFrame, padding="12 12 12 12", relief=RAISED)
        nameFrame.grid(column=0, row=0, columnspan=1, rowspan=16, sticky=(N, W, E, S))
        nameFrame.columnconfigure(0, weight=1)
        nameFrame.rowconfigure(0, weight=1)
        
        # Name entry box
        ttk.Label(nameFrame, text="Name").pack(anchor=W)
        ttk.Entry(nameFrame, textvariable=self.name).pack(anchor=W)

        # Data type radiobuttons
        Radiobutton(master=nameFrame, text="Real data", variable=self.dataType, value=1).pack(anchor=W)
        Radiobutton(master=nameFrame, text="CSEM", variable=self.dataType, value=2).pack(anchor=W)
        Radiobutton(master=nameFrame, text="AxiSEM", variable=self.dataType, value=3).pack(anchor=W)

        ttk.Button(master=nameFrame, text="Load", command=lambda: self.loadData(distCanvas, distAx)).pack(anchor=W)

        # Frequency options and update button
        Radiobutton(master=nameFrame, text="10-30s", variable=self.freqBand, value=1).pack(anchor=W)
        Radiobutton(master=nameFrame, text="10-20s", variable=self.freqBand, value=2).pack(anchor=W)
        Radiobutton(master=nameFrame, text="Custom", variable=self.freqBand, value=3).pack(anchor=W)
        ttk.Label(master=nameFrame, text="Min. Period").pack(anchor=W)
        ttk.Entry(master=nameFrame, textvariable=self.freqMax).pack(anchor=W)
        ttk.Label(master=nameFrame, text="Max. Period").pack(anchor=W)
        ttk.Entry(master=nameFrame, textvariable=self.freqMin).pack(anchor=W)

        ttk.Button(master=nameFrame, text="Update", command=lambda: self.updateFreqBand(distCanvas, distAx)).pack(anchor=W)        

        # Add some padding to all widgets in this frame
        for child in nameFrame.winfo_children():
            child.pack_configure(padx=5, pady=5)

        # Create frame for deconvolution options
        nameFrame2 = ttk.Frame(mainFrame, padding="12 12 12 12", relief=RAISED)
        nameFrame2.grid(column=1, row=0, columnspan=1, rowspan=16, sticky=(N, W, E, S))
        nameFrame2.columnconfigure(0, weight=1)
        nameFrame2.rowconfigure(0, weight=1)

        ttk.Label(nameFrame2, text="Deconvolution options").pack(anchor=W)        

        # Data type radiobuttons
        Radiobutton(master=nameFrame2, text="Water level", variable=self.deconType, value=1).pack(anchor=W)
        Radiobutton(master=nameFrame2, text="Iterative", variable=self.deconType, value=2).pack(anchor=W)

        ttk.Label(master=nameFrame2, text="Water level").pack(anchor=W)
        ttk.Entry(master=nameFrame2, textvariable=self.waterLevel).pack(anchor=W)
        ttk.Label(master=nameFrame2, text="Max. Bumps").pack(anchor=W)
        ttk.Entry(master=nameFrame2, textvariable=self.maxBumps).pack(anchor=W)        
        
        ttk.Button(master=nameFrame2, text="Deconvolve", command=lambda: self.deconvolve()).pack(anchor=W)

        # Add some padding to all widgets in this frame
        for child in nameFrame2.winfo_children():
            child.pack_configure(padx=5, pady=5)        

        
        # Create frame for source plot

        # Create frame for listbox
        listFrame = ttk.Frame(mainFrame, padding="12 12 12 12", relief=RAISED)
        listFrame.grid(column=2, row=0, columnspan=1, rowspan=22, stick=(N, W, E, S))
        listFrame.columnconfigure(0, weight=1)
        listFrame.rowconfigure(0, weight=1)

        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode=SINGLE)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
   
        self.listBox.bind('<Double-Button-1>', self.onDouble)

        # Create canvas for distance plot
        #distFig = plt.figure(figsize=(5,9), dpi=100)
        #self.distAx  = distFig.add_axes([0.1, 0.1, 0.8, 0.8])
        #self.distCanvas = FigureCanvasTkAgg(distFig, master=mainFrame)
        #self.distCanvas.get_tk_widget().grid(row=0, column=9, rowspan=22, columnspan=2)
        #self.distCanvas.show()

        # Create canvas for source plot
        sourceFig = plt.figure(figsize=(10,5), dpi=100)
        self.sourceAx  = sourceFig.add_axes([0.1, 0.1, 0.8, 0.8])
        self.sourceCanvas = FigureCanvasTkAgg(sourceFig, master=mainFrame)
        self.sourceCanvas.get_tk_widget().grid(row=0, column=4, rowspan=16, columnspan=8)
        self.sourceCanvas.show()

        root.bind("<Return>", lambda _: self.loadData())

    def loadData(self):
        dataType = self.dataType.get()
        dir = 'Data/' + str(self.name.get())
        if (dataType == 1):
            self.dir = dir + '/RealData/'
        if (dataType == 2):
            self.dir = dir + '/CSEM/'
        if (dataType == 3):
            self.dir = dir + '/AxiSEM/'

        dirLength = len(self.dir)

        self.seislist = glob.glob(self.dir + '*PICKLE')

        stationList = []
        for x in range(len(self.seislist)):
            station = self.seislist[x]
            station = station[dirLength:]
            station = station[:(len(station) - 7)]
            stationList.append(station)
            
        print(stationList)
        
        self.insertSeisList(stationList)

    def labelUpdate(self):
        self.distance.set(str(round(self.dist, 3)))

    def insertSeisList(self, seislist):
        for item in seislist:
            self.listBox.insert(END, item)

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

        self.plot(canvas, ax)

    def plot(self, canvas, ax):
        # Clear current plot
        ax.clear()

        # Get current data
        seis = read(self.seislist[int(self.index)], format='PICKLE')

        seis.filter('highpass', freq=self.fmin, corners=2, zerophase=True)
        seis.filter('lowpass', freq=self.fmax, corners=2, zerophase=True)

        self.tshift = seis[0].stats['starttime'] - seis[0].stats['eventtime']

        self.seisT = seis.select(channel='BHT')[0]
        
        self.seisT.data = np.gradient(self.seisT.data, self.seisT.stats.delta)

        norm = np.max(abs(self.seisT.data))
        self.seisT.data = self.seisT.data / norm

        self.seisT.times = self.seisT.times() + self.tshift - seis[0].stats.traveltimes['S'] 

        ax.plot(self.seisT.times, self.seisT.data, 'k', linewidth=2)

        ax.set_title('Distance : %.2f' % (seis[0].stats['dist']), loc='left')

        ax.set_ylim([-1.1, 1.1])
        ax.set_xlim([-50, 100])

        if (self.sourceWindow):
            #canvas.create_rectangle(self.x1, 1, self.x2, -1)
            self.sourceWindow = False

        cid = canvas.mpl_connect('button_press_event', self.__onclick__)
            
        canvas.draw()

    def deconvolve(self):
        if (self.idx1 == 0):
            print("You need to select a source first")
        else:
            deconvolve_data(self.dataType.get(),self.name.get(), self.value, self.idx1, self.idx2, self.deconType.get())
        
                    
    def selectSource(self, name, idx1, idx2):

        xArr = self.seisT.times
        array = self.seisT.data

        windowTimes = xArr[idx1:idx2]
        windowData = array[idx1:idx2]

        sourceData = open(self.dir + 'sourceData_'+ name + '.txt', 'w')
        sourceData.write("%s \n" % name)
        sourceData.write("%s \n" % idx1)
        sourceData.write("%s \n" % idx2)
        sourceData.close()

        self.sourceWindow = True

        self.updateFreqBand(self.sourceCanvas, self.sourceAx)


    def onDouble(self, event):
        # Note here that Tkinter passes an event object to onselect()
        widget = event.widget
        self.index = int(widget.curselection()[0])
        self.value = widget.get(self.index)
        self.clickCounter = 0

        self.updateFreqBand(self.sourceCanvas, self.sourceAx)

    def __onclick__(self, click):
        print('Click!')
        self.clickCounter += 1
        self.point = (click.xdata, click.ydata)
        xArr = self.seisT.times
        array = self.seisT.data

        # Find the index of the click in the time array
        if (self.clickCounter == 1):
            self.idx1 = (np.abs(xArr - click.xdata)).argmin()
            self.x1 = click.xdata
            print("Click!", self.idx1, self.x1)

        if (self.clickCounter == 2):
            self.idx2 = (np.abs(xArr - click.xdata)).argmin()
            self.x2 = click.xdata
            print("Click!", self.idx2, self.x2)
            # Window the source
            self.selectSource(self.value, self.idx1, self.idx2)
            
    def instructions(self):
        toplevel = Toplevel()
        label1 = Label(toplevel, text=self.ABOUT_TEXT, height=0, width=100)
        label1.pack(anchor=W)
        label2 = Label(toplevel, text=self.DISCLAIMER, height=0, width=100)
        label2.pack()        

#------------------------------------------------------------------------------------------#
### Initialise class ----------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

root = tk.Tk()
root.title("Part III Project - Source Selection for Deconvolution")
app = Application(master=root)
app.mainloop()

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#    
