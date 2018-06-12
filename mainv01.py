# -*- coding: utf-8 -*-
"""
Created on Wed May 30 12:02:01 2018

@author: Peter
"""

############################################################################################
#                                                                                          #
#   @Morgan okay this is v. messy. Essentially I'm calling                                 #
#   Load                                                                                   #
#   >>wave04, flux04, ln_wave04, ln_flux04, n04 = ld("2004et_20041027_3299_9327_00.dat")   #
#   Process                                                                                #
#   >>filtered_flux04 = pr(ln_wave04,ln_flux04)                                            #
#   correlate                                                                              #
#   >>r = cr(ln_wave04, filtered_flux04,ln_wave04z_half, filtered_flux04)                  # 
#   for every array in order to find how well they correlate to one another                #
#                                                                                          #
############################################################################################
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#from preprocessing import Process as pr
import loader
import preprocessing
import cross_corr
#################################################################################################
ld=loader.Loader()

def Loadin(filename):
    wave,flux = ld.Load(filename)
    bn = loader.Binning(wave)
    ln_flux,ln_wave,N = bn.ln_bin_flux(wave,flux)   
    return wave,flux,ln_wave,ln_flux,N
wave,flux,ln_wave,ln_flux,N = Loadin("2004et_20041027_3299_9327_00.dat")
percent = 0.5
pr = preprocessing.Preproccess(wave,flux,percent)
filtered_sig = pr.Filter()
cr = cross_corr.Correlate(ln_wave,ln_wave,ln_flux,ln_flux)
r_list = cr.correlate()
#wave04, flux04, ln_wave04, ln_flux04, n04 =           ld.Load("2004et_20041027_3299_9327_00.dat")
#wave04p1, flux04p1, ln_wave04p1, ln_flux04p1, n04p1 = ld("d2004et_20041027_3299_9327_00_0p1.dat")
#wave04p2, flux04p2, ln_wave04p2, ln_flux04p2, n04p2 = ld("d2004et_20041027_3299_9327_00_0p2.dat")
#wave04p3, flux04p3, ln_wave04p3, ln_flux04p3, n04p3 = ld("d2004et_20041027_3299_9327_00_0p3.dat")
#wave04z1, flux04z1, ln_wave04z1, ln_flux04z1, n04z1 = ld("2004et_20041027_3299_9327_00_z0p1.dat")

#filtered_flux04 = pr(ln_wave04,ln_flux04)
#filtered_flux04p1 = pr(ln_wave04p1, ln_flux04p1)
#filtered_flux04p2 = pr(ln_wave04p2, ln_flux04p2)
#filtered_flux04p3 = pr(ln_wave04p3, ln_flux04p3)
#filtered_flux04z1 = pr(ln_wave04z1, ln_flux04z1)
#from cross_corr import correlate as cr
#correlation_curve04,
#ln_wave04z_half =[]
#ln_wave04z_half = [element - 0.5 for element in ln_wave04]
#r = cr(ln_wave04, filtered_flux04,ln_wave04z_half, filtered_flux04)
#plt.plot(r)
#correlation_curve04p1, 
#r04p1 = cr(ln_wave04, filtered_flux04,ln_wave04p1, filtered_flux04p1)
###correlation_curve04p2, 
#r04p2 = cr(ln_wave04, filtered_flux04,ln_wave04p2, filtered_flux04p2)
###correlation_curve04p3, 
#r04p3 = cr(ln_wave04, filtered_flux04,ln_wave04p3, filtered_flux04p3)
###correlation_curve04z1, 
#r04z1 = cr(ln_wave04, filtered_flux04,ln_wave04z1, filtered_flux04z1)

#plt.plot(correlation_curve04, label = "04")
#plt.legend()
#plt.figure()
#plt.plot(correlation_curve04p1, label = "04p1")
#plt.legend()
#plt.figure()
#plt.plot(correlation_curve04z1, label = "04z1")
#plt.legend()
#plt.figure()
#plt.plot(correlation_curve04p2, label = "04p2")
#plt.legend()
#plt.figure()
#plt.plot(correlation_curve04p3, label = "04p3")
#plt.legend()
#plt.figure()



######################################################################################
#
#wave11,flux11,ln_wave11,ln_flux11, n11 = ld("sn2011fe_0928_all_cor.txt")
#
#filtered_flux11 = pr(ln_wave11,ln_flux11)
#
##correlation_curve11,
#r11 = cr(ln_wave04, filtered_flux04, ln_wave11, filtered_flux11)
#plt.plot(r11)
#plt.figure()
##plt.plot(correlation_curve11, label = "11")
##plt.legend()
##plt.figure()
#
########################################################################################
##
#wave_99, flux_99, ln_wave_99, ln_flux_99, n_99 =                ld("1999em_19991119_3292_10074_00.dat")
##wave_d99, flux_d99, ln_wave_d99, ln_flux_d99, n_d99 =           ld("d1999em_19991108_4010_10997_00.dat")
##wave_d99p1, flux_d99p1, ln_wave_d99p1, ln_fluxd_99p1, n_d99p1 = ld("d1999em_19991108_4010_10997_00_0p1.dat")
##wave_d99p2, flux_d99p2, ln_wave_d99p2, ln_flux_d99p2, n_d99p2 = ld("d1999em_19991108_4010_10997_00_0p2.dat")
##wave_d99p3, flux_d99p3, ln_wave_d99p3, ln_flux_d99p3, n_d99p3 = ld("d1999em_19991108_4010_10997_00_0p3.dat")
##
#filtered_flux_99 = pr(ln_wave_99,ln_flux_99)
##filtered_flux_99d = pr(ln_wave_d99, ln_flux_d99)
##filtered_flux_99dp1 = pr(ln_wave_d99p1, ln_fluxd_99p1)
##filtered_flux_99dp2 = pr(ln_wave_d99p2, ln_flux_d99p2)
##filtered_flux_99dp3 = pr(ln_wave_d99p3, ln_flux_d99p3)
###
####plt.plot( ln_wave_d99, ln_flux_d99, ln_wave_d99p1, ln_fluxd_99p1, ln_wave_d99p2, ln_flux_d99p2, ln_wave_d99p3, ln_flux_d99p3)
####plt.figure()
###
####correlation_curve99, 
#r99 = cr(ln_wave04, filtered_flux04, ln_wave_99, filtered_flux_99)
#plt.plot(r99)
###correlation_curve99_d, 
#r99_d = cr(ln_wave04, filtered_flux04, ln_wave_d99, filtered_flux_99d)
###correlation_curve99_dp1, 
#r99_dp1 = cr(ln_wave04, filtered_flux04, ln_wave_d99p1, filtered_flux_99dp1)
###correlation_curve99_dp2, 
#r99_dp2 = cr(ln_wave04, filtered_flux04, ln_wave_d99p2, filtered_flux_99dp2)
###correlation_curve99_dp3, 
#r99_dp3 = cr(ln_wave04, filtered_flux04, ln_wave_d99p3, filtered_flux_99dp3)
#
##plt.plot(correlation_curve99, label = "99")
##plt.plot(correlation_curve11, label = "11")
##plt.legend()
##plt.figure()
##plt.plot(ln_wave11,ln_flux11,ln_wave04, ln_flux04)
##plt.figure()
##plt.plot(ln_wave_d99, ln_flux_d99,ln_wave04, ln_flux04)
#
##plt.plot(correlation_curve99_d, label = "99d")
##plt.legend()
##plt.figure()
##plt.plot(correlation_curve99_dp1, label = "99dp1")
##plt.legend()
##plt.figure()
##plt.plot(correlation_curve99_dp2, label = "99dp2")
##plt.legend()
##plt.figure()
##plt.plot(correlation_curve99_dp3, label = "99dp3")
##plt.legend()
##plt.figure()
######################################################################################
#
#wavega,fluxga,ln_wavega,ln_fluxga, n_ga = ld("v_irr_sn_2018.txt")
#ln_wavegaz = ln_wavega+np.log(1+0.4)
##
#filtered_fluxga = pr(ln_wavega,ln_fluxga)
#filtered_fluxgaz = pr(ln_wavegaz,ln_fluxga)
##
###correlation_curvega, 
#rga = cr(ln_wavega, ln_fluxga, ln_wavega, ln_fluxga)
###correlation_curvegaz, 
#rgaz = cr(ln_wavega, ln_fluxga,ln_wavegaz, ln_fluxga)
##
##plt.plot(correlation_curvega, label = "ga")
##plt.legend()
##plt.figure()
##plt.plot(correlation_curvegaz, label = "gaz")
##plt.legend()
##plt.figure()

