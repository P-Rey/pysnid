# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 17:33:21 2018

@author: Peter
"""
import numpy as np

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
