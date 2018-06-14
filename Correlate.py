# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 13:13:40 2018

@author: Peter
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
class Correlate():
    def __init__(self, usr, temp):
        
        self.usrinput = usr
        self.template = temp
        self.start = temp[0]
        self.end = temp[-1]
    
    def get_corr(self):
        inputdft = fft.fft(self.usrinput)
        templatedft = fft.fft(self.template)
#        corr_nondft = inputdft * np.conjugate(templatedft)
#        corr_nondft = np.append(np.zeros(500), corr_nondft)
#        corr = fft.fft(corr_nondft)
        
        self.corr = np.correlate(self.usrinput, self.template, "full")
        CORR = self.corr / ( np.std(inputdft) * np.std(templatedft))
        plt.plot(self.usrinput)
        plt.plot(self.template)
        plt.figure()
        plt.plot(CORR)
        plt.savefig("2004etv2011fe_corr.png")
#        plt.xlim(490,550)
        self.h = max(CORR)
        
        return CORR, self.h
    
    def get_rmsa(self):
        
        auto_corr = np.correlate(self.template, self.template, "full")
        ##### VVVV DIFFERENT LENGTH ISSUE VVVVV #####
        from scipy import interpolate
        
        n = auto_corr[0]
        m = auto_corr[-1]
        x = np.linspace(n, m, len(auto_corr))
        f = interpolate.interp1d(x, auto_corr)

        x1 = np.linspace(n, m, 1024*7)
        f_t = f(x1)
        
        n1 = self.corr[0]
        m1 = self.corr[-1]
        x1 = np.linspace(n1, m1, len(self.corr))
        f1 = interpolate.interp1d(x1, self.corr)
        x2 = np.linspace(n1, m1, 1024*7)
        f_d = f1(x2)
        
        a_n = f_d - f_t
#        a_n = self.corr - auto_corr
        
        rmsa = np.std(a_n)
        
        return rmsa
    
    def get_r(self):
        
        rmsA = self.get_rmsa()
        ### VVVVVV ###
        print("h", self.h)
        print("rmsa", rmsA)
        if rmsA == 0.0:
            r = float(np.inf)
        else:
            r =  self.h / ( np.sqrt(2) * rmsA ) 
        
        return r 
    
    def get_lap(self):
        
        lap = abs ( np.log(self.end/self.start) )
        
        return lap