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

def fitsreadspec():
    spectrum = 1
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

#Plot the spectraBOTTOM
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
filtered_input_spectrum=fftspec(input_spectrum)
filtered_template_spectrum=fftspec(template_spectrum)

# Put onto log wavelength scale
log_input_spectrum = logspec(filtered_input_spectrum)
log_template_spectrum = logspec(filtered_template_spectrum)

# Trim spectra so that we remove unnecessary padding at extrema
log_input_spectrum, log_template_spectrum = logspectrim(log_input_spectrum, log_template_spectrum, input_limits, template_limits)

print("log0",input_spectrum.shape, input_spectrum[0,0])#:,0])#.log_input_spectrum[0])

print(input_spectrum[-1,0])

lap = np.log(log_input_spectrum[-1,0]/log_input_spectrum[0,0])

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
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
## Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
#plt.show()
plt.figure()
#plt.savefig("04et_11fe.pdf")

####################################################################

# Work in progress, find peak in correlation, look at how this compares to the midpoint
maxi =  (np.argmax(correlation))
leng = (len(correlation))

print ((leng/2)-maxi)
print (leng/2)

a_n = np.zeros(len(correlation))

for j in range(1, int(len(correlation))):
    a_n[j] = correlation[-j] - correlation[j-1]

rmsA = np.std(a_n[:int(len(a_n)/2)])

h = max( correlation )
r = h / ( np.sqrt(2) * rmsA )



#lap = np.log(log_input_spectrum[:,-1]/log_input_spectrum[:,0])

rlap =  r * lap

rlap = str(rlap)
rmsA = str(rmsA)
h = str(h)
lap = str(lap)


###################################################################
import tkinter as tk 
import matplotlib.backends.backend_tkagg as tkagg
from tkinter import ttk
from tkinter import *
LARGE_FONT= ("Verdana", 12)
NORM_FONT= ("Verdana", 10)
SMALL_FONT= ("Verdana", 8)


f = Figure()
a = f.add_subplot(111)

def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()
    
def rlappopup():
    root = tk.Tk()
#    popup = tk.Tk()
    root.wm_title("rlaps")
    S = tk.Scrollbar(root)
    T = tk.Text(root, height=4, width=50)
    S.pack()
    T.pack()
    S.config(command=T.yview)
    T.config(yscrollcommand=S.set)
    quote = str("rlap: " + rlap + ", rmsA: " + rmsA + ", h: " + h + ", lap: " + lap)
#    T.insert()
    T.insert(tk.END, quote)
    B1 = ttk.Button(root, text="Okay", command = root.destroy)
    B1.pack()
    root.mainloop()
    
class Application(tk.Tk):

    def __init__(self, *args, **kwargs):
#        window = tk.Tk()
        tk.Tk.__init__(self, *args, **kwargs)

        #tk.Tk.iconbitmap(self, default="clienticon.ico")
        tk.Tk.wm_title(self, "PySNID")
        
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)


        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label = "rlap values", command = lambda: rlappopup())
        filemenu.add_separator()
        filemenu.add_command(label="Save settings", command = lambda: popupmsg("Not supported just yet!"))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="File", menu=filemenu)

        tk.Tk.config(self, menu=menubar)


        fig, ax1  = plt.subplots(1,1)
        ax1.plot(input_spectrum[:,0],input_spectrum[:,1], label = 'User Spectrum')
        ax1.plot(template_spectrum[:,0],template_spectrum[:,1], label = 'Template Spectrum')
        ax1.set_title('Input v Best Fit Template')
        ax1.set_ylabel('Normalised Flux (arb. units)')
        ax1.set_xlabel('Wavelength (Angstrom)')
        ax1.legend()

        
        canvas = FigureCanvasTkAgg(fig, master=container)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=0)
        tkagg.NavigationToolbar2Tk(canvas, container)
        
        OPTIONS = [
        'sn1980K',
        'sn1989B',
        'sn2004et',
        'sn2012aw',
        ] #etc
        
        #master = Tk()
        
        variable = tk.StringVar(container)
        variable.set(OPTIONS[0]) # default value
        
        w = tk.OptionMenu(container, variable, *OPTIONS)
        w.pack()
        def glines(ax1):
            ax1.vlines(x= 6500,ymin =-1,ymax = 2, colors = 'm',linestyles = 'dashdot', label = 'Hydrogen Alpha' )
        b = tk.Button(container, text = 'Galaxy lines', command = lambda:glines(ax1))
        b.pack()
#        self.frames = {}

#        for F in (GraphPage):
#
#            frame = F(container, self)
#
#            self.frames[F] = frame
#
#            frame.grid(row=0, column=0, sticky="nsew")
#
#        self.show_frame(GraphPage)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()

class GraphPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Graph Page!", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(GraphPage))
        button1.pack()

        

        

        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

app = Application()

#app.geometry('1280x720')

app.mainloop()


























####################
'''
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
## Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


###########################

import tkinter as tk
from matplotlib import style
import matplotlib.backends.backend_tkagg as tkagg
LARGE_FONT= ("Verdana", 12)

class PySNID():
    window = tk.Tk()
    window.title("PySNID")
    window.geometry("1270x720")
    
    menubar = tk.Menu()
    filemenu = tk.Menu(menubar)
    filemenu.add_command(label = 'Exit', command = quit)
    
    fig, ax1  = plt.subplots(1,1)
    fig2, ax2 =  plt.subplots(1,1)
    ax1.plot(input_spectrum[:,0],input_spectrum[:,1])
    ax1.plot(template_spectrum[:,0],template_spectrum[:,1])
    ax1.set_title('Input v Best Fit Template')
    ax1.set_ylabel('Normalised Flux (arb. units)')
    ax1.set_xlabel('Wavelength (Angstrom)')
    ax2.plot(correlation)
    
    canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=0)
    tkagg.NavigationToolbar2Tk(canvas, window)
    
    #canvas1 = FigureCanvasTkAgg(fig2, master=window)  # A tk.DrawingArea.
    #canvas1.draw()
    #canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    # canvas is your canvas, and root is your parent (Frame, TopLevel, Tk instance etc.)
    #tkagg.NavigationToolbar2Tk(canvas1, window)
    
    def create_window(window):
        page = tk.Toplevel(window)
    
    def onpress(ax1):
        #if chk_state.set != True:
        #        return
        #else:
        ax1.vlines(x= 6500,ymin =-1,ymax = 2, colors = 'm',linestyles = 'dashdot', label = 'Hydrogen Alpha' )
        #ax1.label()
    chk_state = tk.BooleanVar()
    chk_state.set(True)
    #chk = tk.Checkbutton(window, text ='Galaxy Lines', var = onpress(chk_state, ax1))
    chk = tk.Radiobutton(window, text='Galaxy Lines', variable=onpress(ax1), value=1)
    #chk = tk.Button(window, text ='Galaxy Lines', command = onpress(ax1))
    chk.pack()
    b = tk.Button(window, text= "rlap values", command = create_window(window))
    b.pack()
#    chk.grid(column = 0, row  = 0)
    window.mainloop()
    
app = PySNID()
#app.mainloop()
'''
'''
class SeaofBTCapp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
#        tk.Tk.iconbitmap(self, default ="@/Downloads/ucd_brandmark_colour.ico")
        tk.Tk.wm_title(self, "PySNID")

        container = tk.Frame(self)BOTTOM
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff = 0)
        filemenu.add_command(label = 'Exit', command = quit)
#        filemenu.add_command(label = 'G-Lines', command = onpress())

        self.frames = {}

        frame = StartPage(container, self)
        self.frames[StartPage] = frame
        frame.grid(row=0, column=0, sticky="nsew")
      #  self.show_frame(StartPage)

   # def show_frame(self, cont):

    #    frame = self.frames[cont]
     #   frame.tkraise()


#def qf(quickPrint):
#    print(quickPrint)

        
#class StartPage(tk.Frame):

#    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
#        label = tk.Label(self, text="Start Page", font=LARGE_FONT)
#        label.pack(pady=10,padx=10)
#
#        button = tk.Button(self, text="Visit Page 1",
#                            command=lambda: qf("Check me out, I'm passing vars!"))
#        button.pack()
        root = tk.Tk()
    
        root.wm_title("PySNID")
        
        fig, ax1  = plt.subplots(1,1)
        fig2, ax2 =  plt.subplots(1,1)
        ax1.plot(input_spectrum[:,0],input_spectrum[:,1])
        ax1.plot(template_spectrum[:,0],template_spectrum[:,1])
        ax1.set_title('Input v Best Fit Template')
        ax1.set_ylabel('Normalised Flux (arb. units)')
        ax1.set_xlabel('Wavelength (Angstrom)')
        ax2.plot(correlation)
        
        def onpress():
 #           if event.button != 1:
#                return
            ax1.vlines(x= 6500,ymin =-1,ymax = 2, colors = 'm',linestyles = 'dashdot', label = 'Hydrogen Alpha' )
            ax1.label()
     
        b =tk. Button(root, text = 'Galaxy Lines', command = onpress)
        b.pack()
       # fig.canvas.mpl_connect('button_press_event', onpress)
        
        canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        tkagg.NavigationToolbar2Tk(canvas, root)
        
        canvas1 = FigureCanvasTkAgg(fig2, master=root)  # A tk.DrawingArea.
        canvas1.draw()
        canvas1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
        
        # canvas is your canvas, and root is your parent (Frame, TopLevel, Tk instance etc.)
        tkagg.NavigationToolbar2Tk(canvas1, root)

app = SeaofBTCapp()
app.geometry("1280x720")
#ani = animation.FuncAnimation(fig, animate, interval=1000)

app.mainloop()
'''