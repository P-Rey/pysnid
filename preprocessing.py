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

def Apodize(flux,wave):
    from scipy.interpolate import UnivariateSpline

    p=int(len(wave)/13)                                                      

    b_pntarr = flux[::p]

    a_pntarr = wave[::p] 
    spl = UnivariateSpline(a_pntarr, b_pntarr) 
    spl_b=spl(wave) 
    b_2=flux-spl_b
    Mean=np.mean(b_2)
    SignalSplined=b_2-Mean
    
    return SignalSplined

def Hann(SignalSplined,percent):
    N = len(SignalSplined)
    nsquash = int(percent*N)
    win_len = np.linspace(0,nsquash-1,nsquash-1)
    window = []
    window[:] = [0.5*(1-np.cos(np.pi*element/(nsquash-1))) for element in win_len]
    hanned_sig = SignalSplined[:nsquash-1]*window
    window_rev=window[::-1]
    hanned_sig1 = SignalSplined[-nsquash+1:]*window_rev
    SignalSplined_cut = SignalSplined[:-nsquash+1]
    SignalSplined_cut = SignalSplined_cut[nsquash-1:]
    ProcessedSig = np.append(hanned_sig, SignalSplined_cut)
    ProcessedSig = np.append(ProcessedSig, hanned_sig1)
    
    return ProcessedSig

def Filter(ProcessedSig):
    import scipy.fftpack as fft

    dft=fft.fft(ProcessedSig)
    
    for i in range(50,len(dft)):
        dft[i]=0

    filtered_sig = fft.ifft(dft)

    return filtered_sig

def Process(ln_wave,ln_flux):
    SplinedSignal = Apodize(ln_flux,ln_wave)
    processedsig = Hann(SplinedSignal, 0.05)
    filtered_sig = Filter(processedsig)                                                                  
    return filtered_sig