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
    
        b_pntarr = self.flux[::self.p]
        a_pntarr = self.wave[::self.p] 
        spl = UnivariateSpline(a_pntarr, b_pntarr) 
        spl_b=spl(self.wave) 
        b_2=self.flux-spl_b
        Mean=np.mean(b_2)
        self.SignalSplined=b_2-Mean
        return self.SignalSplined
    
    def Hann(self):
        self.SignalSplined = self.Apodize()
        nsquash = int(self.percent*self.N)
        win_len = np.linspace(0,nsquash-1,nsquash-1)
        window = []
        window[:] = [0.5*(1-np.cos(np.pi*element/(nsquash-1))) for element in win_len]
        hanned_sig = self.SignalSplined[:nsquash-1]*window
        window_rev=window[::-1]
        hanned_sig1 = self.SignalSplined[-nsquash+1:]*window_rev
        SignalSplined_cut = self.SignalSplined[:-nsquash+1]
        SignalSplined_cut = SignalSplined_cut[nsquash-1:]
        self.ProcessedSig = np.append(hanned_sig, SignalSplined_cut)
        self.ProcessedSig = np.append(self.ProcessedSig, hanned_sig1)
        return self.ProcessedSig
    
    def Filter(self):
        self.ProcessedSig = self.Hann()
        dft=fft.fft(self.ProcessedSig)
        
        for i in range(50,len(dft)):
            dft[i]=0
    
        filtered_sig = fft.ifft(dft)
    
        return filtered_sig