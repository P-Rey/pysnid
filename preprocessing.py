# -*- coding: utf-8 -*-
"""
Created on Thu May 17 16:04:41 2018

@author: Peter
"""

################################################################################
#                                                                              #
#   @Morgan this is all for filtering and processing the spectra               #
#   I'm confident its accurate and each function is pretty self explainitory   # 
#   let me know if you need anything clarified                                 #
#                                                                              #
################################################################################

import numpy as np
from scipy.interpolate import UnivariateSpline
import scipy.fftpack as fft

class Preproccess(object):
    
    def __init__(self,wave,flux,percent):
        self.p = int(len(wave)/13)
        self.percent = percent
        self.N = len(wave)
        self.flux = flux
        self.wave = wave
        
    def Apodize(self):
#       apodizes the flux signal using UnivariateSpline 
        b_pntarr = self.flux[::self.p]              
                                                    #   essentially divides each spectra into 13 points     #
        a_pntarr = self.wave[::self.p]              

        spl = UnivariateSpline(a_pntarr, b_pntarr)  #   creates a 13 point spline funciton this is a smooth function that is analogous to the continuum background #
        spl_b=spl(self.wave)                        #   applies this to the wave function                                                                          #
        b_2=self.flux-spl_b                         #   removes the spline (read: continuum background) from the flux data                                         #
        Mean=np.mean(b_2)                           #   finds the mean value of the spline function                                                                #
        self.SignalSplined=b_2-Mean                 #   takes the mean away from the signal, making the mean value zero                                            #
        return self.SignalSplined
    
    def Hann(self):
        self.SignalSplined = self.Apodize()                                             #   gets the apodized signal
        nsquash = int(self.percent*self.N)                                              #   gets the pecetnage of tjhe signal that we want hanned
        win_len = np.linspace(0,nsquash-1,nsquash-1)                                    #   the length of the window to be hanned 
        window = []                                             
        window[:] = [0.5*(1-np.cos(np.pi*element/(nsquash-1))) for element in win_len]  #   creates the Hann function and applies across the length of the desired window
        hanned_sig = self.SignalSplined[:nsquash-1]*window                              #   applies the function to the front half of our signal 
        window_rev=window[::-1]                                                         #   reverses the window function
        hanned_sig1 = self.SignalSplined[-nsquash+1:]*window_rev                        #   applies the reversed window to the back halkf of our signal
        SignalSplined_cut = self.SignalSplined[:-nsquash+1]                             #   
        SignalSplined_cut = SignalSplined_cut[nsquash-1:]                               #   cuts the central region of the signal, which hasn't been hanned
        self.ProcessedSig = np.append(hanned_sig, SignalSplined_cut)
        self.ProcessedSig = np.append(self.ProcessedSig, hanned_sig1)                   #   stitches together all three of the signals. Created a new signal, that is identical to the
                                                                                        #   splined signal, but hanned at smoothed to zeros at either end 
        return self.ProcessedSig
    
    def Filter(self):
        self.ProcessedSig = self.Hann()
        dft=fft.fft(self.ProcessedSig)
        
        for i in range(50,len(dft)):
            dft[i]=0
    
        filtered_sig = fft.ifft(dft)
    
        return filtered_sig