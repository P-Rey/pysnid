# -*- coding: utf-8 -*-
"""
Created on Tue May 22 09:34:07 2018

@author: Peter
"""
import numpy as np
import matplotlib.pyplot as plt
class Correlate(object):
    def __init__(self, data_wave, temp_wave, flux_data, temp_flux):
        ##############################setting limits and initialising arrays###########################
        self.flux_data = flux_data
        self.temp_flux = temp_flux
        self.data_wave = data_wave
        self.temp_wave = temp_wave
        self.N = int(len(self.data_wave))
        print("n",self.N)
        self.Nt = int(len(self.temp_wave))
        self.data_wave_reg = [0] * self.N 
        self.temp_wave_reg = np.zeros(self.Nt+1)
        self.data_flux_reg = np.zeros(self.N)
        self.temp_flux_reg = np.zeros(self.Nt+1)
        self.index_m = np.zeros(self.Nt)
        self.index_n = np.zeros(self.N)
        
    def region_blocked(self,wave_shift):
        
        ###############################################################################################
        #                                                                                             #
        #   this sets limits based on start and end points and it also finds the overlap value lap.   #
        #   The function rebins this into a new array and also cuts out any zeros at either end       #
        #                                                                                             #
        ###############################################################################################

        ##############################assign the start and end points##################################
        
        if wave_shift[0] >= self.data_wave[0]:
            self.start = wave_shift[0]
        else:
            self.start = self.data_wave[0]
        if wave_shift[-1] >= self.data_wave[-1]:
            self.end = self.data_wave[-1]
        else:
            self.end = wave_shift[-1]
            
        self.lap=self.lap_()
       
        ######################################rebinning################################################

        try:
            if self.data_wave[0] > self.start:
                for i in range(0,self.N):
    
                    self.data_wave_reg[i] = self.data_wave[i]
                    self.data_flux_reg[i] = self.flux_data[i]
                    
            
            if wave_shift[0] >=self.start and wave_shift[0] <= self.end:
                print("1",len(self.temp_wave_reg))
                print("2",len(wave_shift))  
                for j in range(0,self.Nt):    
                    self.temp_wave_reg[j] = wave_shift[j]
                    self.temp_flux_reg[j] = self.temp_flux[j]
            
            for k in range(0,self.Nt):
                if self.temp_wave_reg[k] != 0 :
                    self.temp_flux_reg[k] = self.temp_flux[k]
            
            for l in range(0,self.N):
                if self.data_wave_reg[l] != 0:
                    self.data_flux_reg[l] = self.flux_data[l]
        
        #########################deleting zeros at the ends of the new arrays##########################
            for n in range(0,self.N):
                if self.data_wave_reg[n] == 0:
                    self.index_n[n]=n        
                
            for m in range(0,self.Nt):
                if self.temp_wave_reg[m] == 0:
                    self.index_m[m]=m   
            plt.plot(self.index_n)
            plt.figure()
            plt.plot(self.index_m)
#            plt.plot( self.data_wave_reg, self.data_flux_reg, self.temp_wave_reg,self.temp_flux_reg)
#            plt.figure()
            
            self.data_wave_reg = np.delete(self.data_wave_reg, self.index_n)
            self.data_flux_reg = np.delete(self.data_flux_reg, self.index_n)
            
            self.temp_wave_reg = np.delete(self.temp_wave_reg, self.index_m)
            self.temp_flux_reg = np.delete(self.temp_flux_reg, self.index_m)
            
            #plt.plot( self.data_wave_reg, self.data_flux_reg, self.temp_wave_reg,self.temp_flux_reg)
        except IndexError:
            print("AHHHHHHHHH")
        return self.data_wave_reg, self.data_flux_reg, self.temp_wave_reg, self.temp_flux_reg, self.lap
   
    def lap_(self):
        ####################################find lap###################################################
        start_n = np.exp(self.start)
        end_n = np.exp(self.end)
        self.lap=abs(np.log(end_n/start_n))
        return self.lap
    
    def interp_spec(self):
        
        ###################################################################
        #   interpolates any 2 arrays into a new array of length 1024*7   #
        #           because that was found to maximise r                  #
        ###################################################################
        
        from scipy import interpolate
        f = interpolate.interp1d(self.data_wave, self.flux_data)
        n = self.data_wave[0]
        m = self.data_wave[-1]
        x = np.linspace(n, m, 1024*7)
        f_d = f(x)
        
        ft = interpolate.interp1d(self.temp_wave, self.temp_flux)
        nt = self.temp_wave[0]
        mt = self.temp_wave[-1]
        xt = np.linspace(nt, mt, 1024*7)
        f_t = ft(xt)
        return f_d, f_t
    
    def get_r_value(self, h):
        
        #####################################################################################
        #                                                                                   #
        #             a function that takes the spectrum. determines the peak               #
        #           and the rms antisymetrical function of the correlation curve            #
        #                  uses these to find r the correlation coefficient                 #
        #                                                                                   #    
        #####################################################################################
        
        f_d, f_t = self.interp_spec()
        
        
        autocorr1 = np.correlate(f_t , f_t , "full")
        rmsTemp = np.std(autocorr1)
        autocorr1 = (1/len(autocorr1)*rmsTemp*rmsTemp)*autocorr1
    
        Corr1 = np.correlate(f_d, f_t, "full")
        rmsInput = np.std(Corr1)
        Corr1 = (1. / len(Corr1) * rmsInput * rmsTemp) * Corr1
        
        arandom = Corr1- autocorr1
        rmsa = np.std(arandom)
    
        if rmsa == 0:
            import math
            self.r = math.inf
        else:
            self.r = abs((h)/(np.sqrt(2)*rmsa))
        
        return self.r
    
    
    def get_correlate(self, wave_shift):
        
        #############################################################
        #                                                           #
        #   finds the correlation beteern two spectra, fins their   # 
        #   correlation coefficients, calls get_r_value to do so    #
        #   calls region_blocked aswell to find lap                 #
        #                                                           #
        #############################################################
        
        import scipy.fftpack as fft
        
        '''
        SORT OUT THE R COEFFICIENT, STILL TOO LOW. LOOK INTO THE NORMALISATION
        '''
       
        wave_data, flux_data, temp_wave, temp_flux, lap = self.region_blocked(wave_shift)
        plt.figure()
        plt.plot( wave_data, flux_data, temp_wave, temp_flux)
#        dft_data = fft.fft(flux_data)
#        dft_temp = fft.fft(temp_flux)
#        
#        rmsInput = np.std(dft_data)
#        rmsTemp = np.std(dft_temp)
        
        Corr = np.correlate(flux_data, temp_flux, "full")
#        Corr = (1. / len(Corr) * rmsInput * rmsTemp) * Corr
        h = max(Corr)
        r = self.get_r_value(h)
        rlap=r*lap
    
#            rlap = 0
            
        return rlap
    
    def correlate(self):
        
        ################################################################################################    
        #                                                                                              #
        #   This function loops through many redshifts, calls rbl from redshift.py and get_correlate   #
        #                                                                                              #
        ################################################################################################
    
        from redshift import red_boi_loop as rbl
        
        r_list = []
          
        for i in range(0,300):
            wave_shift = rbl(self.temp_wave,i)
            r = self.get_correlate(wave_shift)
            r_list = np.append(r_list, r)
    
        return r_list