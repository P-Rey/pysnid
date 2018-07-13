# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:28:32 2018

@author: Peter
"""

import numpy as np
import pandas as pd

class Loader(object):
    def Load(self, filename):
        
        #########################################
        #                                       #
        #   loads and Normalises  the spectra   #
        #                                       #
        #########################################
        
        '''
        Clean this up and make it more consistent. Also why am I using DataFrame here, what is it good for?
        '''
        data = np.loadtxt(filename)
        data1 = pd.DataFrame(data)
        wave = data1.loc[:, 0]
        flux = data1.loc[:, 1]
        wave=np.asarray(wave)
        flux = np.asarray(flux)

        return wave, flux
    def template_loader(self,direct, filename):
        template_file = np.genfromtxt(os.path.join(direct, filename))
class Binning(object):
    def __init__(self, wave):
        
        ##############################################
        #                                            #
        #   Gets the constants from the wave array   #
        #                                            #
        ##############################################
        
        self.l0=wave[0]
        self.l1=wave[-1]
        self.N=len(wave)
        
        self.l10 = self.l1/self.l0
        self.A = self.N / np.log(self.l10)
        
        self.B = -self.N * (np.log(self.l0)/np.log(self.l10))
    
        self.dl_ln = np.log(self.l10) / self.N 
        self.ln_wave = np.log(wave)
    
        self.n= []
        self.n[:] = [self.A * np.log(element) + self.B for element in wave]
    
    
    def ln_bin_flux(self, wave, flux):
        
        ####################################################
        #                                                  #
        #                 From apodize.f                   #  
        #                                                  #
        ####################################################
    
        '''
        can this be vecotrised?
        '''
        ln_flux = np.zeros(self.N)
        for i in range(0,self.N):
            if i == 0:
                s0 = 0.5 * (3 * wave[i] - wave[i + 1])
                s1 = 0.5 * (wave[i] + wave[i + 1])
            elif i == len(wave) - 1:
                s0 = 0.5 * (wave[i - 1] + wave[i])
                s1 = 0.5 * (3 * wave[i] - wave[i - 1])
            else:
                s0 = 0.5 * (wave[i - 1] + wave[i])
                s1 = 0.5 * (wave[i] + wave[i + 1])
                
            ln_s0= np.log(s0 / self.l0) / self.dl_ln + 1
            ln_s1 = np.log(s1 / self.l0) / self.dl_ln + 1
            dnu = s1 - s0
            
            for j in range(int(ln_s0), int(ln_s1)):
                if j < 0 or j >= self.N:
                    continue
                flux_ = flux[i] / (ln_s1 - ln_s0) * dnu
                ln_flux[j] = ln_flux[j] + flux_
#        ln_flux = (ln_flux - min(ln_flux)) / (max(ln_flux) - min(ln_flux))
        return ln_flux, self.ln_wave, self.N



#################################################################################
'''
#ln_flux = np.zeros(int(N)) #empty array to be filled with the final ln binned flux 
#ln_wave = np.log(wave)
#ln_bin = np.zeros(int(N))
#for i in range(1,N):
#    ln_bin[i-1]= ln_wave[i]-ln_wave[i-1]
#maybe_lin_bin = np.exp(ln_bin)
##ln_bin_widths = maybe_lin_bin
##for i in range(0,N):
##    for j in range(0,N): 
##        if  np.log(wave[i])<= ln_bin[j+1] & np.log(wave[i])>=ln_bin[j]      :
##            some_thing_happens
#hbs = np.zeros(N)
#lbs = np.zeros(N)
#ln_hbs = np.zeros(N)
#ln_lbs = np.zeros(N)
#
#for i in range(1,N):      # for flux[i] 
#    hbs[i-1] = wave[i]-0.5*dl_ln # higher_bin_side
#    lbs[i] = wave[i-1]+0.5*dl_ln # lower_bin_side
#    
#    ln_hbs[i-1] = ln_wave[i]-0.5*maybe_lin_bin[i] 
#    ln_lbs[i] = ln_wave[i-1]+0.5*maybe_lin_bin[i-1]  
#
#hbs = np.log(hbs[:-1])
#lbs = np.log(lbs[:1])

hbs can be >=<ln_hbs and >=< ln_lbs
lbs can be >=<ln_hbs and >=< ln_lbs
but hbs > lbs and ln_hbs > ln_lbs

find the where hbs[i] is closest ln_hbs and lbs[i] is closest ln_lbs[i]
is hbs >=< ln_hbs? 

if hbs > lnhbs then lbs has all options 

if hbs = ln_hbs then lbs !>= ln_hbs

if hbs < ln_hbs then lbs !>=< ln_hbs

elif hbs>ln_lbs but <ln_hbs lbs =< ln_lbs # maybe check ln_lbs[i-1]
def lbs_check(newflux[i],lbs,ln_lbs):
    if lbs[i] == ln_lbs[i]:
        newflux[i] = newflux[i]
    if lbs[i] > ln_lbs[i]:
        newflux[i] = newflux[i]+flux[i-1]
    if lbs[i] < ln_lbs[i]:
        #flux[i] is contributing to flux[i-1] this would have been picked up in an earlier loop
        
#####Assumption made: seeing as the spacing is never significantly different(hbs isnt >> or << ln_hbs) I assum that hbs is always > ln_hbs
#newflux = np.zeros(N)
#for j in range(0,N-1):
#    if j == 0:
#        continue
#    elif hbs[j] >= ln_hbs[j]:                    #this means that flux[i+1] is also in this bin
#        newflux[j] = flux[j]             
#        if lbs[j] > ln_lbs[j]:
#            newflux[j] = 0.5*(newflux[j]+flux[j-1])
#    elif hbs[j] < ln_hbs[j]:                    #this means that flux[i] carrys into the next bin
#        if j == N-1:
#            newflux[i] = flux[i]
#        else:
#            newflux[j] = 0.5(flux[j]+flux[j+1])                    #assumes that 50% of both fluxes contribute.this doesn't account for if the hbs[i] < lbs[i]
#            if lbs[j] > ln_lbs[j]:
#                newflux[j] = 0.5*(newflux[j]+flux[j-1])    


#maybe_lin_bin[0] = maybe_lin_bin[0]+l0
#maybe_ = np.zeros(N)
#for j in range(1,N):
#    maybe_lin_bin[j] = maybe_lin_bin[j] + maybe_lin_bin[j-1]
#maybe_ = maybe_lin_bin+l0

#if wave[i]<maybe_[i]<wave[i+1]
#waveMiddle = wave[1:-1] #everything but the first and last wavelengths             
                          
#np.log(some_list/ l0) / dl_ln + 1


first bin starts at =l0
las bin ends at at l1
second bin is at 

nonln wave bin 
ln_wavebin

#############
'''
