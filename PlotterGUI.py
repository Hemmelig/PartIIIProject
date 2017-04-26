#!/usr/bin/env python3

import sys
import os
import subprocess
from tkinter import *
import tkinter as tk
from tkinter import ttk
import obspy
from obspy import read, UTCDateTime
from obspy.core import Stream
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
sys.path.insert(0, '/raid2/cb793/GUI')
from plot_data_with_distance import *

class Application(tk.Frame):
    ABOUT_TEXT = """Instructions for use:

    This software can be used to download and process data from the IRIS
    catalogue website. """
    
    DISCLAIMER = """Author

    This interface was developed by Me."""

    # Filter frequencies
    fmin = 0.033
    fmax = 0.1
    
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        
        self.name = StringVar()
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
        
        # Frequency options and update button
        Radiobutton(master=optionFrame, text="10-30s", variable=self.freqBand, value=1).pack(anchor=W)
        Radiobutton(master=optionFrame, text="10-20s", variable=self.freqBand, value=2).pack(anchor=W)
        Radiobutton(master=optionFrame, text="Custom", variable=self.freqBand, value=3).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Min. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMax).pack(anchor=W)
        ttk.Label(master=optionFrame, text="Max. Period").pack(anchor=W)
        ttk.Entry(master=optionFrame, textvariable=self.freqMin).pack(anchor=W)

        # Control buttons for plots     
        ttk.Button(master=optionFrame, text="Data: Distance", command=lambda i=1: self.plot(i)).pack(anchor=W)
        ttk.Button(master=optionFrame, text="Data: Azimuth", command=lambda i=2: self.plot(i)).pack(anchor=W)
        ttk.Button(master=optionFrame, text="CSEM: Azimuth", command=lambda i=3: self.plot(i)).pack(anchor=W)        
        ttk.Button(master=optionFrame, text="Map: Turning Points", command=lambda i=4: self.plot(i)).pack(anchor=W)        
        ttk.Button(master=optionFrame, text="Map: Delay Times", command=lambda i=5: self.plot(i)).pack(anchor=W)        
        ttk.Button(master=optionFrame, text="2D Delay Time Array", command=lambda i=6: self.plot(i)).pack(anchor=W)
        
        # Add some padding to all widgets in this frame
        for child in optionFrame.winfo_children():
            child.pack_configure(padx=5, pady=5)        

        self.update()

    def plot(self, i):
        # Clear current plot
        ax.clear()
        n = str(self.name.get())
        print(i, n)

        if (i == 1):
            plotDataDistance(n, self.fmin, self.fmax)
        elif (i == 2):
            plotDataAzimuth(n, self.fmin, self.fmax)
        elif (i == 3):
            plotCSEMAzimuth(n, self.fmin, self.fmax)
        elif (i == 4):
            plotMapTurn(n, self.fmin, self.fmax)
        elif (i == 5):
            plotMapDelay(n, self.fmin, self.fmax)
        elif (i == 6):
            plot2Ddt(n, self.fmin, self.fmax)
        
        canvas.draw()

    def instructions(self):
        toplevel = Toplevel()
        label1 = Label(toplevel, text=self.ABOUT_TEXT, height=0, width=100)
        label1.pack(anchor=W)
        label2 = Label(toplevel, text=self.DISCLAIMER, height=0, width=100)
        label2.pack()        
    
root = tk.Tk()
root.title("Part III Project - Plot Creator")
app = Application(master=root)
app.mainloop()
