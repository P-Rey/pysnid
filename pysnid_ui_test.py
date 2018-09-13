#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 11:23:52 2018

@author: peter
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate, fftpack, signal
import seaborn as sns
import matplotlib.animation as animation

#matplotlib inline
#plt.rcParams['figure.figsize'] = [14, 8]
#sns.set()
#sns.set(rc={'figure.figsize':(14, 8)})

# Function to read in spectrum from ascii file
# Returns 2d numpy with wavelength, flux columns
# TO DO - Check if ascii file is valid, strip out header
#       - Write equivalent function for fits files

def readspec(filename):
    # Read in spectrum
    spectrum = []
    with open(filename) as file:
        for line in file:
            x = line.split()
            spectrum.append([float(x[i]) for i in range(len(x))])
    # Convert to numpy array
    spectrum = np.array(spectrum)
    return spectrum

def parsereadspec():
    # Read in spectrum
    parser = argparse.ArgumentParser()
    parser.add_argument("input_spec")
    args = parser.parse_args()
    spectrum = []
    with open(args.input_spec) as file:
        for line in file:
            x = line.split()
            spectrum.append([float(x[i]) for i in range(len(x))])
    # Convert to numpy array
    spectrum = np.array(spectrum)
    return spectrum

# Function to resample spectrum onto grid of log wavelengths
# We have a grid from ln(lambda) = 7.0 to 10.0, with 0.001 steps
# This has 3000 bins between approx 1100 Ang to 2.2 um, covering 
# the full UVOIR region. We have 1000 bins between ~3000 to 9000
# Ang, which is well matched to the typical sampling of low 
# resolution SN spectra.
#
# The function returns a 2d numpy array containing:
# [ln(lambda), lambda, flux]

def logspec(spectrum):
    # Set up empty array
    logspectrum = np.zeros([3000,3])

    # Create a 1d array with the log wavelength values
    logwave = np.arange(len(logspectrum[:,1]))
    logwave = 7+logwave/1000

    # Copy this into the first colum of our output array
    logspectrum[:, 0] = logwave
    # And copy the equivalent 
    logspectrum[:, 1] = np.exp(logspectrum[:,0])


    logspectrum[:,2] = np.interp(logspectrum[:,1], spectrum[:,0], spectrum[:,1], left=0, right=0)
    
    return logspectrum

# Function to Fourier filter spectra
def fftspec(spectrum):

    fft = fftpack.fft(spectrum[:,1])
    filtered_fft = fft.copy()
    filtered_fft[100:]=0    # Can play around with this value (100) to filter more or less
    filtered_spectrum = fftpack.ifft(filtered_fft)
    spectrum[:,1] = filtered_spectrum
    
    return spectrum

# Function to fit a spline to a spectrum and return a spectrum normalised
# about zero
def splinespec(spectrum):

    spline_fit = interpolate.UnivariateSpline(spectrum[:,0],spectrum[:,1])
    normalised_spectrum = (spectrum[:,1]/spline_fit(spectrum[:,0]))-1

    spectrum[:,1] =  normalised_spectrum
    return spectrum

# Function to return shortest and longest wavelengths in spectrum
# Output of this is required by logspectrim
def getspeclimits(spectrum):
    spec_lambda_blue = spectrum[0,0]
    spec_lambda_red  = spectrum[-1,0]

    return [spec_lambda_blue, spec_lambda_red]


# Function to trim log spectra
# At red end we just trim to the shortest spectrum.
# At blue end we trim to shortest wavelength (i.e. we allow for the input spectrum to move to the blue)
def logspectrim(spec1, spec2, spec1limits, spec2limits):

    if spec1limits[1] == spec2limits[1]:
        lambda_max = spec1limits[1]
    else:
        lambda_max = min(spec1limits[1], spec2limits[1])

    # Since our log wavelength grid is at 0.001 steps, we convert the wavelength
    # to a log wavelength with this precision
    ln_lambda_max = np.around(np.log(lambda_max), decimals=3)
    # Now search for the index of the spectrum which corresponds to this wavelength
    index_lambda_max = np.where(spec1[:,0]==ln_lambda_max)[0][0]

    # Need to find index_lambda_min, index of min (spec1[0,1], spec2[0,1])
    if spec1limits[0] == spec2limits[0]:
        lambda_min = spec1limits[0]
    else:
        lambda_min = min(spec1limits[0], spec2limits[0])
    ln_lambda_min = np.around(np.log(lambda_min), decimals=3)
    index_lambda_min = np.where(spec1[:,0]==ln_lambda_min)[0][0]

    spec1 = spec1[index_lambda_min:index_lambda_max,:]
    spec2 = spec2[index_lambda_min:index_lambda_max,:]

    
    return spec1, spec2

# Custom filter
# work in progress, maybe better with a cosine tail at the ends

def mfilt (spectrum):
    filter = np.ones(len(spectrum))
    length =len(spectrum)
    for i in range(0,int(length*0.1)):
        filter[i] = i/(length*0.1)
    for i in range(int(length*0.90),length):
        filter[i] =  (float(length-i)/((length-(length*0.90))+1))
        
    return filter    

# Read in files for template and input spectrum
import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument("input")
#args = parser.parse_args()
input_spectrum = parsereadspec()



#template_spectrum=readspec('Desktop/2004et_20041027_3299_9327_00.dat')
template_spectrum = readspec('SN2011fe_2011-08-25_00-00-00_TNG_DOLORES_PTF.ascii')
#input_spectrum=readspec('2004et_20041027_3299_9327_00.dat')

#Plot the spectra
plt.plot(input_spectrum[:,0], input_spectrum[:,1])
plt.plot(template_spectrum[:,0], template_spectrum[:,1])

#plt.show()
plt.figure()
##########################################

# Find longest and shortest wavelengths in each spectrum
input_limits = getspeclimits(input_spectrum)
template_limits = getspeclimits(template_spectrum)

# Normalise spectra with spline fit
input_spectrum=splinespec(input_spectrum)
template_spectrum=splinespec(template_spectrum)

# Fourier filter out high frequency noise
input_spectrum=fftspec(input_spectrum)
template_spectrum=fftspec(template_spectrum)

# Put onto log wavelength scale
log_input_spectrum = logspec(input_spectrum)
log_template_spectrum = logspec(template_spectrum)

# Trim spectra so that we remove unnecessary padding at extrema
log_input_spectrum, log_template_spectrum = logspectrim(log_input_spectrum, log_template_spectrum, input_limits, template_limits)


plt.plot(log_input_spectrum[:,0], log_input_spectrum[:,2])
plt.plot(log_template_spectrum[:,0], log_template_spectrum[:,2])

#############################################################

# Make a Hann filter
hann = signal.hann(len(log_template_spectrum))
#filter = mfilt(log_template_spectrum)

# Apply Hann filter to spectra
log_template_spectrum = log_template_spectrum[:,2] * hann
log_input_spectrum = log_input_spectrum[:,2] * hann

plt.plot(log_template_spectrum)
plt.plot(log_input_spectrum)
plt.figure()
correlation = signal.correlate(log_input_spectrum, log_template_spectrum, mode='full')

plt.plot(correlation, label="")
plt.legend()
#plt.show()
plt.figure()
#plt.savefig("04et_11fe.pdf")

####################################################################

# Work in progress, find peak in correlation, look at how this compares to the midpoint
maxi =  (np.argmax(correlation))
leng = (len(correlation))

print ((leng/2)-maxi)
print (leng/2)

###################################################################

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
## Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
#from six.moves import tkinter as Tk
#
#root = Tk.Tk()
#root.wm_title("Embedding in Tk")
#
#fig = Figure(figsize=(5, 2), dpi=100)
#t = np.arange(0, 3, .01)
#fig.add_subplot(111).plot( correlation)
#fig.add_subplot(111).plot( log_template_spectrum)
#canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
#canvas.draw()
#canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
#
#toolbar = NavigationToolbar2TkAgg(canvas, root)
#toolbar.update()
#canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
#
#
#def on_key_press(event):
#    print("you pressed {}".format(event.key))
#    key_press_handler(event, canvas, toolbar)
#
#
#canvas.mpl_connect("key_press_event", on_key_press)
#
#
#def _quit():
#    root.quit()     # stops mainloop
#    root.destroy()  # this is necessary on Windows to prevent
#                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#
#
#button = Tk.Button(master=root, text="Quit", command=_quit)
#button.pack(side=Tk.BOTTOM)
#
#Tk.mainloop()

###########################

import tkinter as tk
from matplotlib import style
import matplotlib.backends.backend_tkagg as tkagg
LARGE_FONT= ("Verdana", 12)
#style.use("ggplot")
#
#fig = Figure(figsize=(10, 4), dpi=100)
#a = fig.add_subplot(111)
#
#
#def animate(i):
#    pullData = open('sampleText.txt','r').read()
#    dataArray = pullData.split('\n')
#    xar=[]
#    yar=[]
#    for eachLine in dataArray:
#        if len(eachLine)>1:
#            x,y = eachLine.split(',')
#            xar.append(int(x))
#            yar.append(int(y))
#    a.clear()
#    a.plot(xar,yar)
    
class SeaofBTCapp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
#        tk.Tk.iconbitmap(self, default ="@/Downloads/ucd_brandmark_colour.ico")
        tk.Tk.wm_title(self, "PySNID")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}



        frame = StartPage(container, self)

        self.frames[StartPage] = StartPage(container, self) # frame

        frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()


def qf(quickPrint):
    print(quickPrint)

        
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
#        label = tk.Label(self, text="Start Page", font=LARGE_FONT)
#        label.pack(pady=10,padx=10)
#
#        button = tk.Button(self, text="Visit Page 1",
#                            command=lambda: qf("Check me out, I'm passing vars!"))
#        button.pack()
        root = tk.Tk()
    
       # root.wm_title("PySNID")
        fig = Figure(figsize=(10, 4), dpi=100)
        fig1 = Figure(figsize=(10, 4), dpi=100)

        fig.add_subplot(111).plot(input_spectrum[:,0],input_spectrum[:,1])
        
        fig.add_subplot(111).plot(template_spectrum[:,0],template_spectrum[:,1])
        fig1.add_subplot(121).plot(correlation)
        
#        figsrc, axsrc = plt.subplots()
#        figzoom, axzoom = plt.subplots()
#        axsrc.set(xlim=(0, 1), ylim=(0, 1), autoscale_on=False,title='Click to zoom')
#        axzoom.set(xlim=(0.45, 0.55), ylim=(0.4, 0.6), autoscale_on=False,title='Zoom window')
#
#        def onpress(event):
#            if event.button != 1:
#                return
#            x, y = event.xdata, event.ydata
#            axzoom.set_xlim(x - 0.1, x + 0.1)
#            axzoom.set_ylim(y - 0.1, y + 0.1)
#            figzoom.canvas.draw()
#        
#        figsrc.canvas.mpl_connect('button_press_event', onpress)
        
        canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        tkagg.NavigationToolbar2Tk(canvas, root)
        
        canvas1 = FigureCanvasTkAgg(fig1, master=root)  # A tk.DrawingArea.
        canvas1.draw()
        canvas1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        
        # canvas is your canvas, and root is your parent (Frame, TopLevel, Tk instance etc.)
        tkagg.NavigationToolbar2Tk(canvas1, root)

app = SeaofBTCapp()

#ani = animation.FuncAnimation(fig, animate, interval=1000)

app.mainloop()