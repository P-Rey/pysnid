# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:53:06 2018

@author: Peter
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
#from itertools import izip

text = np.loadtxt(os.path.join("templates-2.0_2500_20000_1024", "sn2004et.lnw"))
#file = pd.read_csv(os.path.join("templates-2.0_2500_20000_1024", "sn2004et.lnw"), sep='\s+' )
#read_csv  sn98dx_bsnip  

#for i, sub_list in enumerate(file):
#    for j, value in enumerate(sub_list):
        #print('a[{}][{}] = {}'.format(i, j, value))
#print(file[1:])
wave=text[1:,0]
flx=text[1:,30]
N = len(text)
listoboy = np.linspace(0,N, N)
z = text[:1]
inter = iter(z)
#from itertools import izip
b = dict(zip(listoboy,z))
plt.plot(wave,flx)
plt.xlim(2700,7700)
#a=[]
#for i in range(len(file)):
#    if i == float('NaN'):
#        a[i] = i
#print(next(i for i in file if i == 0.0))
    
#file = open(os.path.join("bsnip_v7_snid_templates.tar", "sn04et_bsnip.lnw"), "r")
#mytxt = file.readlines()
#file.close()
#
#del mytxt[0]
#n = len(mytxt[0])
#N = len(mytxt)

#new_txt = []
#for i in range(0,N):
#    if len(mytxt[i]) != n:
#        new_txt = np.append(new_txt, mytxt[i])
##NN = len(new_txt)
#
##np.split(new_txt," ")
##print(new_txt[1])
##value = dict()
#value = "value"
#for x in range(1,len(new_txt)+1):
#    value+str(x) = np.zeros(len(new_txt))
    #my_array = np.zeros(len(new_txt))
#n_df = pd.Series(new_txt)
#n_N_df = n_df.str.split()
#NNN = len(n_N_df)
#NN = len(n_N_df[0])
#for m in range(NN):
#    str("array"+m) = np.zeros(NN)
#for j in range(NNN): 
#    for k in range(NN):
#        n_N_df[j][k]

#for i in range(20):
#    print(len(new_dataframe[i]))
#print(FAVUHB)
#new_file = np.genfromtxt(os.path.join("templates-2.0_2500_20000_2048", "sn2004et.lnw"), skip_header = 12 )
#a =0.3
#wave = new_file[: ,0]
#flux1 = new_file[:,1]
#flux2 = new_file[:,2]+a
#flux3 = new_file[:,3]+a*2
#flux4 = new_file[:,4]+a*3
#flux5 = new_file[:,5]+a*4
#flux6 = new_file[:,6]+a*5
#flux7 = new_file[:,7]+a*6
#flux8 = new_file[:,8]+a*7
#flux9 = new_file[:,9]+a*8
#flux10 = new_file[:,10]+a*9
#flux11 = new_file[:,11]+a*10
#
#plt.plot(wave[1:],flux1[1:])
#plt.plot(wave[1:],flux2[1:])
#plt.plot(wave[1:],flux3[1:])
#plt.plot(wave[1:],flux4[1:])
#plt.plot(wave[1:],flux5[1:])
#plt.plot(wave[1:],flux6[1:])
#plt.plot(wave[1:],flux7[1:])
#plt.plot(wave[1:],flux8[1:])
#plt.plot(wave[1:],flux9[1:])
#plt.plot(wave[1:],flux10[1:])
#plt.plot(wave[1:],flux11[1:])
#plt.xlim(3000,9500)