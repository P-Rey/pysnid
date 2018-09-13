# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 17:34:35 2018

@author: Peter
"""
from scipy import interpolate, fftpack, signal
import numpy as np
# Function to Fourier filter spectra

#class apodiser():
#    
##    def __init__(self, spectrum):
##        self.spectrum = spectrum
        
def fftspec(spectrum):

    fft = fftpack.fft(spectrum[:,1])
    filtered_fft = fft.copy()
    filtered_fft[100:]=0    # Can play around with this value (100) to filter more or less
    filtered_spectrum = fftpack.ifft(filtered_fft)
    spectrum[:,1] = filtered_spectrum
    
    return spectrum

# Function to fit a spline to a spectrum and return a spectrum normalised
# about zero
def splinespec(apodised_spectrum):

    spline_fit = interpolate.UnivariateSpline(apodised_spectrum[:,0],apodised_spectrum[:,1])
    normalised_spectrum = (apodised_spectrum[:,1]/spline_fit(apodised_spectrum[:,0]))-1

    apodised_spectrum[:,1] =  normalised_spectrum
    return apodised_spectrum

def mfilt(apodised_spectrum):
    filter = np.ones(len(apodised_spectrum))
    length =len(apodised_spectrum)
    for i in range(0,int(length*0.1)):
        filter[i] = i/(length*0.1)
    for i in range(int(length*0.90),length):
        filter[i] =  (float(length-i)/((length-(length*0.90))+1))
        
    return filter    

def hanwin(apodised_spectrum):
    endlength = int(0.05*len(apodised_spectrum))
    cosine = signal.hann(endlength*2)
    hann = np.ones(len(apodised_spectrum))
    hann[0:endlength] = cosine[0:endlength]
    hann[-endlength:] = cosine[endlength:endlength*2]
    hann = apodised_spectrum * hann
    return hann

def apodise(spectrum):
    apodised_spectrum = splinespec(spectrum)
    apodised_spectrum = fftspec(apodised_spectrum)
    apodised_spectrum = mfilt(apodised_spectrum)
    apodised_spectrum = hanwin(apodised_spectrum)
    return apodised_spectrum

