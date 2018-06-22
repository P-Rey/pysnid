# -*- coding: utf-8 -*-
"""
Created on Mon May 28 11:33:55 2018

@author: Peter
"""
import numpy as np
from preprocessing import Hann

def DFT(ProcessedSig):
    import scipy.fftpack as fft

    dft=fft.fft(ProcessedSig)
    
    for i in range(50,len(dft)):
        dft[i]=0

    sigma_squared= np.sqrt((1/len(dft))*np.sum((dft)**2))

    return dft, sigma_squared

def Correlate(dftdata,dfttemp,sigma,sigmatemp):
    dftdata=dftdata
    dfttemp=dfttemp
    Nd = int(len(dftdata))

    Init_Correlate = dftdata[:int(len(dfttemp))]*np.conj(dfttemp)
    
    import scipy.fftpack as fft
    Trans_Corr=fft.ifft(Init_Correlate )
    percent= 0.50
    hanned_corr = Hann(Trans_Corr,percent)
    rmsInput = np.std(dftdata)
    rmsTemp = np.std(dfttemp)
    hanned_corr= (1. / Nd * rmsInput * rmsTemp) * hanned_corr   
    return hanned_corr