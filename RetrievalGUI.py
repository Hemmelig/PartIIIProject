#!/usr/bin/env python3

### The purpose of this code is to provide a simple GUI allowing the user to perform
### a range of different actions on a given dataset
### To use:
### 1. Fill in the details of the event and stations of interest
### 2. Download the data then process it
### 3. Add the CMT source data - open base file, edit then save in relevant Data folder
### 4. Add synthetics and traveltimes
### 5. Plot data. !!! Need to add distance and azimuth functionality of plots !!!

#####################################################################################
### Import modules and external python scripts
#####################################################################################

import sys, os, subprocess
from tkinter import *
import tkinter as tk
from tkinter import ttk, filedialog
from tkinter.scrolledtext import ScrolledText
sys.path.insert(0, '/raid2/cb793/GUI')
from download_data_one_event import *
from data_processing import *
from add_synthetics_topickle import *
from add_travel_times import *
from add_event_cmt import *
from plot_data_with_azimuth import *
from plot_data_with_distance import *
from add_turn_points import *

#####################################################################################
    
class Application(tk.Frame):
    ABOUT_TEXT = """Instructions for use:

    This software can be used to download and process data from the IRIS
    catalogue website. """
    
    DISCLAIMER = """Author

    This interface was developed by Me."""

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.name = StringVar()
        self.latitude = StringVar()
        self.longitude = StringVar()
        self.starttime = StringVar()
        self.endtime = StringVar()
        self.maxrad = StringVar()
        self.minmag = StringVar()
        self.maxmag = StringVar()
        self.distmin = StringVar()
        self.distmax = StringVar()
        self.azmin = StringVar()
        self.azmax = StringVar()
        self.lengthoftrace = StringVar()

        self.createWidgets()

    def createWidgets(self):
        # Create menu
        menubar = Menu(root)
        menubar.add_command(label="Instructions", command=lambda: self.instructions())
        menubar.add_command(label="Exit", command=root.quit)
        root.config(menu=menubar)

        # Create frames and populate with widgets
        # Mainframe
        self.mainFrame = ttk.Frame(root, padding="3 3 12 12")
        self.mainFrame.grid(column=0, row=0, sticky=(N, W, E, S))
        self.mainFrame.columnconfigure(0, weight=1)
        self.mainFrame.rowconfigure(0, weight=1)

        # Name and processing button frame
        nameFrame = ttk.Frame(self.mainFrame, padding="12 12 12 12", relief=RAISED)
        nameFrame.grid(column=0, row=0, columnspan=2, rowspan=8, sticky=(N, W, E, S))
        nameFrame.columnconfigure(0, weight=1)
        nameFrame.rowconfigure(0, weight=1)

        ttk.Label(nameFrame, text="Name").pack(anchor=W)
        ttk.Entry(nameFrame, textvariable=self.name).pack(anchor=W)
        ttk.Button(nameFrame, text="1. Download", command=lambda: self.download()).pack(anchor=W)
        ttk.Button(nameFrame, text="2. Process", command=lambda: self.process()).pack(anchor=W)
        ttk.Button(nameFrame, text="3. Add CMT", command=lambda: self.addCMT()).pack(anchor=W)
        ttk.Button(nameFrame, text="4. Add travel times", command=lambda: self.addTravelTimes()).pack(anchor=W)
        ttk.Button(nameFrame, text="5. Add synthetics", command=lambda: self.addSynths()).pack(anchor=W)
        ttk.Button(nameFrame, text="6. Add turn points", command=lambda: self.addTurnPoints()).pack(anchor=W)

        for child in nameFrame.winfo_children():
            child.pack_configure(padx=5, pady=5)         

        # CMT frame
        cmtFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        cmtFrame.grid(column=2, row=0, columnspan=8, rowspan=8, sticky=(N, W, E, S))
        cmtFrame.columnconfigure(0, weight=1)
        cmtFrame.rowconfigure(0, weight=1)

        cmtEntryFrame = ttk.Frame(cmtFrame, padding="3 3 12 12")
        cmtEntryFrame.pack(side=LEFT)
    
        ttk.Label(cmtEntryFrame, text="CMT Source Data").pack(anchor=N)
        
        self.cmtEntry = ScrolledText(cmtEntryFrame)
        self.cmtEntry.pack(anchor=W)

        cmtSaveFrame = ttk.Frame(cmtFrame, padding="12 12 12 12")
        cmtSaveFrame.pack(side=LEFT)

        ttk.Button(cmtSaveFrame, text="Open", command=lambda: self.openCommand()).pack(anchor=W)
        ttk.Button(cmtSaveFrame, text="Save", command=lambda: self.saveCommand()).pack(anchor=W)        

        # Event parameter frame
        eventFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        eventFrame.grid(column=0, row=8, sticky=(N, W, E, S))
        eventFrame.columnconfigure(0, weight=1)
        eventFrame.rowconfigure(0, weight=1)

        ttk.Label(eventFrame, text="Event Parameters").pack(anchor=W)

        # Label names
        labels = ["Latitude", "Longitude", "Start Time", "End Time", "Max Radius", "Min Magnitude", "Max Magnitude"]
        parameters = [self.latitude, self.longitude, self.starttime, self.endtime, self.maxrad, self.minmag, self.maxmag]

        xtmp = zip(labels, parameters)

        for x in xtmp:
            ttk.Label(eventFrame, text=x[0]).pack(anchor=W)
            ttk.Entry(eventFrame, textvariable=x[1]).pack(anchor=W)
            
        del xtmp
        del labels
        del parameters

        # Station parameter frame
        stationFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        stationFrame.grid(column=2, row=8, sticky=(N, W, E, S))
        stationFrame.columnconfigure(0, weight=1)
        stationFrame.rowconfigure(0, weight=1)

        ttk.Label(stationFrame, text="Station Parameters").pack(anchor=W)

        # Label names
        labels = ["Min Distance", "Max Distance", "Min Azimuth", "Max Azimuth", "Length of trace"]
        parameters = [self.distmin, self.distmax, self.azmin, self.azmax, self.lengthoftrace]

        xtmp = zip(labels, parameters)

        for x in xtmp:
            ttk.Label(stationFrame, text=x[0]).pack(anchor=W)
            ttk.Entry(stationFrame, textvariable=x[1]).pack(anchor=W)
            
        del xtmp            

        # Travel time frame
        travelFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        travelFrame.grid(column=4, row=8, sticky=(N, W, E, S))
        travelFrame.columnconfigure(0, weight=1)
        travelFrame.rowconfigure(0, weight=1)

        # Create label, checkbuttons and execute button
        ttk.Label(travelFrame, text="Travel times").pack(anchor=W)
        ttk.Label(travelFrame, text="Select phases:").pack(anchor=W)

        # Loop through phase options and create check buttons for each. Just add more names
        # to the phases list to produce more buttons. May need to add a counter and if
        # statement in order to split into two columns beyond a certain number of phases
        self.phases = ["S", "Sdiff", "ScS", "SKS", "SKKS"]
        self.vPhases = [None] * len(self.phases)

        for i in range(len(self.phases)):
            self.vPhases[i] = BooleanVar()
            Checkbutton(travelFrame, text=self.phases[i], variable=self.vPhases[i]).pack(anchor=W)     
        
        # Turning points frame
        turnFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        turnFrame.grid(column=6, row=8, sticky=(N, W, E, S))
        turnFrame.columnconfigure(0, weight=1)
        turnFrame.rowconfigure(0, weight=1)

        # Create label, checkbuttons and execute button
        ttk.Label(turnFrame, text="Turning points").pack(anchor=W)
        ttk.Label(turnFrame, text="Select phases:").pack(anchor=W)

        # Loop through phase options and create check buttons for each. Just add more names
        # to the phases list to produce more buttons. May need to add a counter and if
        # statement in order to split into two columns beyond a certain number of phases
        self.phases2 = ["S", "ScS", "Sdiff"]
        self.vPhases2 = [None] * len(self.phases2)

        for i in range(len(self.phases2)):
            self.vPhases2[i] = BooleanVar()
            Checkbutton(turnFrame, text=self.phases2[i], variable=self.vPhases2[i]).pack(anchor=W)     

        # Plotting frame
        plotFrame = ttk.Frame(self.mainFrame, padding="3 3 12 12", relief=RAISED)
        plotFrame.grid(column=8, row=8, sticky=(N, W, E, S))
        plotFrame.columnconfigure(0, weight=1)
        plotFrame.rowconfigure(0, weight=1)

    def download(self):
        # Get event parameters
        n = str(self.name.get())
        l1 = float(self.latitude.get())
        l2 = float(self.longitude.get())
        st = str(self.starttime.get())
        et = str(self.endtime.get())
        mr = float(self.maxrad.get())
        mim = float(self.minmag.get())
        mam = float(self.maxmag.get())

        # Get station parameters
        dmin = float(self.distmin.get())
        dmax = float(self.distmax.get())
        amin = float(self.azmin.get())
        amax = float(self.azmax.get())
        ltrace = float(self.lengthoftrace.get())

        # Pass arguments to download data script
        downloadData(n, l1, l2, st, et, mr, mim, mam, dmin, dmax, amin, amax, ltrace)

    def process(self):
        # Get event name
        n = str(self.name.get())
        processData(n)   

    def addSynths(self):
        # Get event name
        n = str(self.name.get())
        addSynthetics(n, True)
          
    def read(self, phase):
        return self.vPhases[phase].get()

    def addTravelTimes(self):
        # Get event name
        n = str(self.name.get())

        # Check for any phases to include
        phasesToInclude = []
        for i in range(len(self.vPhases)):
            if self.read(i):
                phasesToInclude.append(self.phases[i])
            else:
                continue

        travelTimes(n, phasesToInclude)

    def read2(self, phase):
        return self.vPhases2[phase].get()

    def addTurnPoints(self):
        # Get event name
        n = str(self.name.get())

        # Check for any phases to include
        phasesToInclude2 = []
        for i in range(len(self.vPhases2)):
            if self.read2(i):
                phasesToInclude2.append(self.phases2[i])
            else:
                continue

        turningPoints(n, phasesToInclude2)

    def openCommand(self):
        myfile = filedialog.askopenfile(parent=self.mainFrame,
                              title='Open a Python file', mode='r')
        loadedfile = myfile.read()
        myfile.close()
        self.cmtEntry.delete('1.0', END)
        self.cmtEntry.insert("end", loadedfile)

    def saveCommand(self):
        file = filedialog.asksaveasfile(mode='w')
        if file != None:
        # Slice off the last character from get, as an extra return is added
            data = self.cmtEntry.get('1.0', END+'-1c')
            file.write(data)
            file.close()

    def addCMT(self):
        # Get event name
        n = str(self.name.get())   

        addEventCMT(n)

    def instructions(self):
        toplevel = Toplevel()
        label1 = Label(toplevel, text=self.ABOUT_TEXT, height=0, width=100)
        label1.pack(anchor=W)
        label2 = Label(toplevel, text=self.DISCLAIMER, height=0, width=100)
        label2.pack()

#####################################################################################
### Plotter section
#####################################################################################

### Want to extend this section to give better choices: amplitude, distance divisions
### for azimuth plotter etc

#    def plotter(self):
#        # Get event name
#        n = str(name.get())
 #       if (plotType.get() == 1):
 #           plotWithAzimuth(n, syn)
 #       elif (plotType.get() == 2):
 #           plotWithDistance(n, syn)

#plotType = IntVar()
#syn = IntVar()
#Radiobutton(mainframe, text="Azimuth", variable=plotType, value=1).grid(column=4, row=17, columnspan=2, sticky=W)
#Radiobutton(mainframe, text="Distance", variable=plotType, value=2).grid(column=4, row=18, columnspan=2, sticky=W)
#Radiobutton(mainframe, text="Yes", variable=syn, value=1)
#Radiobutton(mainframe, text="No", variable=syn, value=2)

#ttk.Button(mainframe, text="Plot", command=plotter).grid(column=5, row=18, columnspan=2, sticky=W)

#####################################################################################
        
root = tk.Tk()
root.title("Part III Project - Data Retrieval Interface")
app = Application(master=root)
app.mainloop()
