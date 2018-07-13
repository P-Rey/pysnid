# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 15:53:06 2018

@author: Peter
"""
import numpy as np
import os
import yaml
import Correlate

class template():
    
    def __init__(self, wave_usr, flux_usr):
        self.filepath = "template_directory.yml"
        self.directory_location = "templates-2.0_2500_20000_1024"
        self.wave_usr = wave_usr
        self.flux_usr = flux_usr
        self.sn_array = []
        
    def template_loader(self):

        with open(self.filepath, "r") as file_descriptor:
            
            data = yaml.load(file_descriptor)
            
        for i in range(1, len(data)):

            path = "supernova_"+str(i)
            items = data.get(path)
            junk, rhs= items.split("label:'", 1)
            sn_type, junk = rhs.split("type:", 1)
            sn_type, junk = sn_type.split("'", 1)
            del junk
            self.sn_array = np.append(self.sn_array,sn_type)
#            print('here?')
        return self.sn_array
    
    def template_findr(self):
#        print('what?')
        array_names = self.template_loader()
#        print(len(array_names))
        for i in range(len(array_names)):
            
            text = np.loadtxt(os.path.join(str(self.directory_location), str(array_names[i])+'.lnw'))
#            print(os.path.join(str(self.directory_location), str(array_names[i])+'.lnw'))
            text = text.T
            wave = text[0]
            rlap_array = np.zeros(len(text))
            rlap_top = np.zeros(len(array_names))
            print(len(text)) #new for each template
            for j in range(1,len(text)):
                
                flux = text[j]
                corr = Correlate.Correlate(self.wave_usr, self.flux_usr, wave,flux)
                correlate, h, lap = corr.get_corr()
                r = corr.get_r()
                rlap = r * lap
                print(rlap)
                rlap_array = np.append(rlap_array, rlap)
                return rlap_array
            rlap_top = np.append(rlap_top, rlap_array)
            print(len(rlap_top))
            return rlap_top
        return rlap_top


#p=1
#
#flux_1 = text[1] + p * 1
#flux_2 = text[2] + p * 2
#flux_3 = text[3] + p * 3
#flux_4 = text[4] + p * 4
#flux_5 = text[5] + p * 5
#flux_6 = text[6] + p * 6
#flux_7 = text[7] + p * 7
#flux_8 = text[8] + p * 8
#flux_9 = text[9] + p * 9
#flux_10 = text[10] + p * 10
#flux_11 = text[11] + p * 11
#flux_12 = text[12] + p * 12
#flux_13 = text[13] + p * 13
#flux_14 = text[14] + p * 14
#flux_15 = text[15] + p * 15
#flux_16 = text[16] + p * 16
#flux_17 = text[17] + p * 17
#flux_18 = text[18] + p * 18
#flux_19 = text[19] + p * 19
#flux_20 = text[20] + p * 20
#flux_21 = text[21] + p * 21
#flux_22 = text[22] + p * 22
#flux_23 = text[23] + p * 23
#flux_24 = text[24] + p * 24
#flux_25 = text[25] + p * 25
#flux_26 = text[26] + p * 26
#flux_27 = text[27] + p * 27
#flux_28 = text[28] + p * 28
#flux_29 = text[29] + p * 29
#flux_30 = text[30] + p * 30 
#flux_31 = text[31] + p * 31
#flux_32 = text[32] + p * 32
#flux_33 = text[33] + p * 33
#flux_34 = text[34] + p * 34
#flux_35 = text[35] + p * 35
#flux_36 = text[36] + p * 36
#flux_37 = text[37] + p * 37
#flux_38 = text[38] + p * 38
#flux_39 = text[39] + p * 40
#flux_40 = text[40] + p * 41
#flux_41 = text[41] + p * 42
#flux_42 = text[42] + p * 43
#flux_43 = text[43] + p * 44
#flux_44 = text[44] + p * 45
#flux_45 = text[45] + p * 46
#flux_46 = text[46] + p * 47
#flux_47 = text[47] + p * 48
#flux_48 = text[48] + p * 49

#plt.plot(wave[1:], flux_1[1:])
#plt.plot(wave[1:], flux_2[1:])
#plt.plot(wave[1:], flux_3[1:])
#plt.plot(wave[1:], flux_4[1:])
#plt.plot(wave[1:], flux_5[1:])
#plt.plot(wave[1:], flux_6[1:])
#plt.plot(wave[1:], flux_7[1:])
#plt.plot(wave[1:], flux_8[1:])
#plt.plot(wave[1:], flux_9[1:])
#plt.plot(wave[1:], flux_10[1:])
#plt.plot(wave[1:], flux_11[1:])
#plt.plot(wave[1:], flux_12[1:])
#plt.plot(wave[1:], flux_13[1:])
#plt.plot(wave[1:], flux_14[1:])
#plt.plot(wave[1:], flux_15[1:])
#plt.plot(wave[1:], flux_16[1:])
#plt.plot(wave[1:], flux_17[1:])
#plt.plot(wave[1:], flux_18[1:])
#plt.plot(wave[1:], flux_19[1:])
#plt.plot(wave[1:], flux_20[1:])
#plt.plot(wave[1:], flux_21[1:])
#plt.plot(wave[1:], flux_22[1:])
#plt.plot(wave[1:], flux_23[1:])
#plt.plot(wave[1:], flux_24[1:])
#plt.plot(wave[1:], flux_25[1:])
#plt.plot(wave[1:], flux_26[1:])
#plt.plot(wave[1:], flux_27[1:])
#plt.plot(wave[1:], flux_28[1:])
#plt.plot(wave[1:], flux_29[1:])
#plt.plot(wave[1:], flux_30[1:])
#plt.plot(wave[1:], flux_31[1:])
#plt.plot(wave[1:], flux_32[1:])
#plt.plot(wave[1:], flux_33[1:])
#plt.plot(wave[1:], flux_34[1:])
#plt.plot(wave[1:], flux_35[1:])
#plt.plot(wave[1:], flux_36[1:])
#plt.plot(wave[1:], flux_37[1:])
#plt.plot(wave[1:], flux_38[1:])
#plt.plot(wave[1:], flux_39[1:])
#plt.plot(wave[1:], flux_40[1:])
#plt.plot(wave[1:], flux_41[1:])
#plt.plot(wave[1:], flux_42[1:])
#plt.plot(wave[1:], flux_43[1:])
#plt.plot(wave[1:], flux_44[1:])
#plt.plot(wave[1:], flux_45[1:])
#plt.plot(wave[1:], flux_46[1:])
#plt.plot(wave[1:], flux_47[1:])
#plt.plot(wave[1:], flux_48[1:])
#plt.xlim(3500,9000)
#plt.plot(wav[1:]e, flux_)
#[1:]
#flux_49 = text[49]


#fluxtest = np.reshape(flux[:1],(49,))

#for i in range(1,N):
#    plt.plot(wave[-1:],flux[-1:, i])
#    plt.figure()

#wave=text[1:,0]
#flx=text[:,30]
#N = len(text)
#
#
#listoboy = np.linspace(0,N, N)
#listoboy = [int(element) for element in listoboy]
#text_trans = text.T
#flux = text_trans[1:]
#z = text[:1]
#z = z.T
#z = np.reshape(z,(49,))
#inter = iter(z)
#
#zippidy = zip(listoboy,z)
#b = dict(zippidy)
#plt.plot(wave,flx)
#plt.xlim(2700,10000)
#print(b[1])
#print('flux_'+str(b[1]))



#################################################################################
#plt.savefig('snid_temp_04et.png')
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