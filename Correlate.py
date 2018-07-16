# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 13:13:40 2018

@author: Peter
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
from scipy import interpolate

class Correlate():
    def __init__(self, usrwave, usr, templatewave, temp):
        self.usrinput = usr
        self.usrwave = usrwave
        self.template = temp
        self.template_wave = templatewave
        self.start = temp[0]
        self.end = temp[-1]
    
    def region(self):
        
        if self.usrwave[0] > self.template_wave[0]:
            idx_start_template = min(range(len(self.template_wave)), key = lambda i: abs(self.template_wave[i] - self.usrwave[0]))
            idx_start_usr = 0
        else:
            idx_start_usr =  min(range(len(self.usrwave)), key = lambda i: abs(self.usrwave[i] - self.template_wave[0]))
            idx_start_template = 0
            
        if self.usrwave[-1] > self.template_wave[-1]:
            idx_end_usr = min(range(len(self.usrwave)), key = lambda i: abs(self.usrwave[i] - self.template_wave[-1]))
            idx_end_template = -1
        else:
            idx_end_template =  min(range(len(self.template_wave)), key = lambda i: abs(self.template_wave[i] - self.usrwave[-1]))
            idx_end_usr = -1
            
        new_flux = self.usrinput[idx_start_usr:idx_end_usr]
        new_wave = self.usrwave[idx_start_usr:idx_end_usr]

        new_flux_template = self.template[idx_start_template:idx_end_template]
        new_wave_template = self.template_wave[idx_start_template:idx_end_template]
        lap = self.get_lap(new_wave)

        return new_flux, new_flux_template, lap
    
    def get_corr(self):
        new_flux,new_flux_template, lap = self.region()
        inputdft = fft.fft(new_flux)
        templatedft = fft.fft(new_flux_template)
        
        self.corr = np.correlate(new_flux, new_flux_template)
        self.CORR = self.corr / ( np.std(inputdft) * np.std(templatedft))

        self.CORR = np.real(self.CORR)
        self.h = max(self.CORR)
        
        return self.CORR, self.h, lap
    
    def get_rmsa(self):
        auto_corr = np.correlate(self.template,self.template, "full")
        auto_corr = auto_corr/(np.std(fft.fft(self.template)) * np.std(fft.fft(self.template)))

        np.seterr(divide = 'ignore', invalid = 'ignore')
        np.warnings.filterwarnings('ignore')
        n = auto_corr[0]
        m = auto_corr[-1]
        x = np.linspace(n, m, len(auto_corr))
        f = interpolate.interp1d(x, auto_corr)

        x1 = np.linspace(n, m, 1024 * 7)
        f_t = f(x1)
        
        n1 = self.corr[0]
        m1 = self.corr[-1]
        x1 = np.linspace(n1, m1, len(self.CORR))
        f1 = interpolate.interp1d(x1, self.CORR)
        x2 = np.linspace(n1, m1, 1024 * 7)
        f_d = f1(x2)
        a_n = f_d - f_t
        
        autocorr = np.zeros(len(self.CORR))
        
        for i in range(1, len(self.CORR)):
            autocorr[i-1] = self.CORR[-i] - self.CORR[i-1]
        autocorr = np.real(autocorr)
        rmsa = np.std(autocorr)

        return rmsa
    
    def get_r(self):
        
        rmsA = self.get_rmsa()

        if rmsA == 0.0:
            r = float(np.inf)
        else:
            r =  self.h / ( np.sqrt(2) * rmsA ) 
        
        return r 
    
    def get_lap(self, wave):
        
        lap = abs ( np.log(wave[-1]/wave[0]) )
        
        return lap