# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 15:07:49 2018

@author: Peter
"""
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
#data = np.loadtxt("sn1979C.lnw")
op = open("sn1981B.lnw")
listyboi = op.readlines()
N = 2057-8
wave = np.zeros(0)
flux = np.zeros(0)
flux1 = np.zeros(0)
flux2 = np.zeros(0)
flux3 = np.zeros(0)
flux4 = np.zeros(0)
for i in range(11,2056):
    gd = listyboi[i]
    ahhhh= gd.split( )
    wave = np.append(wave,ahhhh[0])
    flux = np.append(flux,ahhhh[1])
    flux1 = np.append(flux1,ahhhh[2])
    flux2 = np.append(flux2,ahhhh[3])
    flux3 = np.append(flux3,ahhhh[4])
    flux4 = np.append(flux4,ahhhh[5])
wave = wave.astype(np.float)
flux = flux.astype(np.float)
flux1 = flux1.astype(np.float)
flux2 = flux2.astype(np.float)
flux3 = flux3.astype(np.float)
flux4 = flux4.astype(np.float)
plt.plot(wave,flux, label = "10")
plt.figure()
plt.plot(wave,flux1,label = "11")
plt.figure()
plt.plot(wave,flux2,label = "13")
plt.figure()
plt.plot(wave,flux3,label = "43")
plt.figure()
plt.plot(wave,flux4, label = "72")
plt.figure()
plt.legend()
plt.ylim(-0.5,1)
plt.xlim(2000,8000)
#data1 = pd.DataFrame("sn1979C.lnw")
#wave = data1.loc[:, 0]
#flux = data1.loc[:, 1]
#wave=np.asarray(wave)
#flux = np.asarray(flux)
#
