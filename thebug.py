import numpy as np
from scipy.interpolate import UnivariateSpline
import scipy.fftpack as fft
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
from math import sqrt

class Loader(object):
    def Load(filename):

        data = np.loadtxt(filename)
        data1 = pd.DataFrame(data)
        wave = data1.loc[:, 0]
        flux = data1.loc[:, 1]
        wave=np.asarray(wave)
        flux = np.asarray(flux)

        return wave, flux
        
class Binning(object):
    def __init__(self, wave, flux):
        self.flux = flux
        self.wave = wave
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
    
    
    def ln_bin_flux(self):
        
        ln_flux = np.zeros(self.N)
        for i in range(0,self.N):
            if i == 0:
                s0 = 0.5 * (3 * self.wave[i] - self.wave[i + 1])
                s1 = 0.5 * (self.wave[i] + self.wave[i + 1])
            elif i == len(self.wave) - 1:
                s0 = 0.5 * (self.wave[i - 1] + self.wave[i])
                s1 = 0.5 * (3 * self.wave[i] - self.wave[i - 1])
            else:
                s0 = 0.5 * (self.wave[i - 1] + self.wave[i])
                s1 = 0.5 * (self.wave[i] + self.wave[i + 1])
                
            ln_s0= np.log(s0 / self.l0) / self.dl_ln + 1
            ln_s1 = np.log(s1 / self.l0) / self.dl_ln + 1
            dnu = s1 - s0
            
            for j in range(int(ln_s0), int(ln_s1)):
                if j < 0 or j >= self.N:
                    continue
                flux_ = self.flux[i] / (ln_s1 - ln_s0) * dnu
                ln_flux[j] = ln_flux[j] + flux_
        ln_flux = (ln_flux - min(ln_flux)) / (max(ln_flux) - min(ln_flux))
        return self.ln_wave, ln_flux, self.dl_ln

class Preproccess(object):
    
    def __init__(self,wave,flux,percent):
        self.p = int(len(wave) / 13 )#/3.1)
        self.percent = percent
        self.N = len(wave)
        self.flux = flux
        self.wave = wave
        
    def Apodize(self):
        
        self.flux = self.flux - np.mean(self.flux)
        b_pntarr = self.flux[::self.p]              
                                                    
        a_pntarr = self.wave[::self.p]              
        print(a_pntarr)
        spl = UnivariateSpline(a_pntarr, b_pntarr)  
        spl_b=spl(self.wave)                        
        b_2=self.flux-spl_b                         
        self.SignalSplined = b_2
        self.SignalSplined = self.SignalSplined - np.mean(self.SignalSplined)
        
        return self.SignalSplined
    
    def Hann(self):
        
        self.SignalSplined = self.Apodize()                                         
        nsquash = int(self.percent*self.N)                                             
        win_len = np.linspace(0,nsquash-1,nsquash-1)                                    
        window = []                                             
        window[:] = [0.5*(1-np.cos(np.pi*element/(nsquash-1))) for element in win_len]  
        hanned_sig = self.SignalSplined[:nsquash-1]*window                            
        window_rev=window[::-1]                                                        
        hanned_sig1 = self.SignalSplined[-nsquash+1:]*window_rev                        
        SignalSplined_cut = self.SignalSplined[:-nsquash+1]                             
        SignalSplined_cut = SignalSplined_cut[nsquash-1:]                              
        self.ProcessedSig = np.append(hanned_sig, SignalSplined_cut)
        self.ProcessedSig = np.append(self.ProcessedSig, hanned_sig1)          
        
        return self.ProcessedSig
    
    def Filter(self):
        self.ProcessedSig = self.Hann()
        dft=fft.fft(self.ProcessedSig)

        for i in range(50,len(dft)):
            dft[i]=0
    
        filtered_sig = fft.ifft(dft)
        return filtered_sig

def Loadin(filename):
    wave,flux = Loader.Load(filename)
    wave = wave/1.0001
    bn = Binning(wave, flux)
    ln_flux,ln_wave,dwlog= bn.ln_bin_flux()   
    return wave,flux,ln_wave,ln_flux, dwlog
#person = input('Enter your name: ')
#sn_file = input('enter filename:! ')
wavep01,fluxp01,ln_wavep01,ln_fluxp01, dwlog = Loadin("2004et_20041027_3299_9327_00.dat")
'sn1994I-19940405.txt'
'1996al_19960910_fina.txt'
'sn2007kj-20071008.txt'
"2004et_20041027_3299_9327_00.dat"
'sn2011fe_0928_all_cor.txt'
'sn12aw_20120517_BC_300_1_1DF'
'sn2004fu-20041110.17-fast'
'1999dn_19990828.txt'
'1990B_19900123_4200_8197_00.txt'
def Preprommm(wave,flux):
    percent = 0.15
    pr = Preproccess(wave,flux,percent)
    filtered_sig = pr.Filter()
    return filtered_sig

filtered_sigp01 = Preprommm(wavep01, fluxp01)
filtered_sigp01= ( filtered_sigp01 - np.mean(filtered_sigp01) ) / ( max(filtered_sigp01) - np.mean (filtered_sigp01) )



'sn1984A'

text = np.loadtxt(os.path.join("templates-2.0_2500_20000_1024", 'sn2004et.lnw')) # 'sn2004et.lnw'
text = text.T
snid_wave = text[0]
snid_wave = snid_wave[1:]
snid_flux = text[13]
snid_flux = snid_flux[1:]
snid_flux = ( snid_flux - np.mean(snid_flux) ) / ( max(snid_flux) - np.mean (snid_flux) )


plt.plot(snid_wave[136:648],snid_flux[136:648])
plt.plot(wavep01,filtered_sigp01)
plt.figure()

def binboy(wave,flux):
    bn = Binning(wave, flux)
    ln_wave,ln_flux,dwlog = bn.ln_bin_flux()
#    filtered_sigp01= ( filtered_sigp01 - np.mean(filtered_sigp01) ) / ( max(filtered_sigp01) - np.mean (filtered_sigp01) )
    ln_flux = (ln_flux - np.mean(ln_flux)) / (max(ln_flux) - np.mean(ln_flux))
    return ln_wave,ln_flux


log_snid_wave,log_snid_flux = binboy(snid_wave, snid_flux)

log_wave,log_flux = binboy(wavep01,filtered_sigp01)

plt.plot(log_snid_wave[136:648],log_snid_flux[136:648])
plt.plot(log_wave,log_flux)
plt.figure()



f = interpolate.interp1d(log_wave,log_flux)
f_snid = interpolate.interp1d(log_snid_wave[136:648],log_snid_flux[136:648])
num = len(log_wave)
x = np.linspace(log_wave[0],log_wave[-1],num)
x_snid = np.linspace(log_snid_wave[136], log_snid_wave[647], num)

#lap = abs()

f = interpolate.interp1d(wavep01, filtered_sigp01)
f_snid = interpolate.interp1d(snid_wave[136:648],snid_flux[136:648])

num = len(wavep01)
x = np.linspace(wavep01[0],wavep01[-1], num)
x_snid = np.linspace(snid_wave[136], snid_wave[647], num)


lap = abs ( np.log(wavep01[-1] / wavep01[0]))

print('lap', lap)

flux_new = f(x)
snid_flux_new = f_snid(x_snid)
plt.figure()

plt.plot(x,flux_new, label = 'input') 
plt.plot(x_snid, snid_flux_new, label ='template')
plt.xlabel('wavelength (Angstrom)')
plt.ylabel('Normalised flux (arb.)')
plt.legend()
plt.savefig('templatevinput.png')
plt.figure()

correlation = np.correlate(flux_new,snid_flux_new, 'full')

correlation = correlation / (np.std(fft.fft(snid_flux)) * np.std(fft.fft(snid_flux_new)) )

autocorr = np.correlate(snid_flux_new, snid_flux_new, 'full')
autcorr = autocorr / (np.std(fft.fft(snid_flux_new)) * np.std(fft.fft(snid_flux_new)) )

a_n = np.zeros(len(correlation))

#correlation = abs(correlation)


for j in range(1, int(len(autocorr))):
    a_n[j] = correlation[-j] - correlation[j-1]
#a_n = correlation - autcorr
rmsA = np.std(a_n[:int(len(a_n)/2)])
middleIndex = (len(correlation) - 1)/2
index_peak = np.argmax(correlation) 


plt.figure()
plt.plot(correlation)
plt.plot(a_n)
#plt.vlines(middleIndex,-0.5,1.5, color = 'r')
#plt.vlines(index_peak,-0.5,1.5)
#plt.xlim(1000,1040)
plt.figure()
autocorr = correlation - a_n

plt.plot(autocorr)
plt.plot(correlation)
print('rmsA', rmsA)
print('argmax',np.argmax(correlation))
shift = abs((2047 /2 )- np.argmax(correlation))
print('shift', shift)
h = max( correlation )
r = h / ( np.sqrt(2) * rmsA )
rlap =  r * lap
#index_peak = np.argwhere(h)

wave_bin = (x[-1] - x[0]) / num
middleIndex = (len(correlation) - 1)/2
bin_shift = middleIndex - index_peak
delta = wave_bin * bin_shift
#dwlog = np.log(x[-1]/x[0]) / len(x)

z = (np.exp( delta*dwlog )-1)
#z = (lamba_obs -lamba_act) / lambda_act 



def calc_redshift_from_crosscorr(crossCorr, nw, dwlog):
    deltaPeak = np.argmax(crossCorr)
    print("deltaPeak",deltaPeak)
    # z = np.exp(deltaPeak * dwlog) - 1 #equation 13 of Blondin)
    zAxisIndex = np.concatenate((np.arange(-nw / 2, 0), np.arange(0, nw / 2)))
    print("zAxisIndex ",zAxisIndex )
    plt.figure()
    plt.plot(zAxisIndex)
    print("nw",nw/2)

    if deltaPeak < nw / 2:
        z = (np.exp(abs(zAxisIndex) * dwlog) - 1)[deltaPeak]
    else:
        z = -(np.exp(abs(zAxisIndex) * dwlog) - 1)[deltaPeak]

    return z
#dwlog = np.log(wavep01[-1]/wavep01[0]) / len(wavep01)
redshift = calc_redshift_from_crosscorr(correlation, len(log_wave), dwlog)

'''
trim_list = np.argwhere(snid_flux == snid_flux[0])
for i in range(1,len(trim_list)):
    if trim_list[i] - trim_list[i-1] != 1:
        trim_high = trim_list[i:]
        trim_low = trim_list[:i]

snid_wave_new = snid_wave[int(trim_low[-1]):int(trim_high[0]-1)]
snid_flux_new = list(filter(lambda x: x != snid_flux[0], snid_flux))

snid_flux_new = np.trim_zeros(snid_flux)

for a in snid_flux_new:
    if a in snid_flux:
#        index_new = np.argwhere( snid_flux_new == a)
        index_snid = np.argwhere( snid_flux == a)
        break
snid_wave_new = snid_wave[int(index_snid[0]) :  ]
smid_wave_new = snid_wave_new[:-int(len(snid_wave) - index_snid[-1] - 1)]#int(snid_wave[index_snid[-1]]):]
######################################################################################################
    
def region(template_wave,template, usrwave, usrinput):

        if usrwave[0] > template_wave[0]:
            idx_start_template = min(range(len(template_wave)), key = lambda i: abs(template_wave[i] - usrwave[0]))
            idx_start_usr = 0
        else:
            idx_start_usr =  min(range(len(usrwave)), key = lambda i: abs(usrwave[i] - template_wave[0]))
            idx_start_template = 0
            
        if usrwave[-1] > template_wave[-1]:
            idx_end_usr = min(range(len(usrwave)), key = lambda i: abs(usrwave[i] - template_wave[-1]))
            idx_end_template = -1
        else:
            idx_end_template =  min(range(len(template_wave)), key = lambda i: abs(template_wave[i] - usrwave[-1]))
            idx_end_usr = -1
            
        new_flux = usrinput[idx_start_usr:idx_end_usr]
        new_wave = usrwave[idx_start_usr:idx_end_usr]

        new_flux_template = template[idx_start_template:idx_end_template]
        new_wave_template = template_wave[idx_start_template:idx_end_template]
        
        if new_wave_template[0] < new_wave[0]:
            new_wave[0] = new_wave_template[0]
        else:
            new_wave_template[0] = new_wave[0]
        if new_wave_template[-1] > new_wave[-1]:
            new_wave[-1] = new_wave_template[-1]
        else:
            new_wave_template[-1] = new_wave[-1]
            
        return new_wave, new_flux, new_wave_template, new_flux_template
    
#new_wave, new_flux, new_wave_template, new_flux_template = region(smid_wave_new,snid_flux_new, wavep01, filtered_sigp01)


######################################################################################################
'''
'''

#############################################################
#                                                           #
#                   Input v TemplatePlot                    #
#                                                           #
#############################################################
plt.figure()
plt.plot(x_snid,snid_flux_new, label = '(template)')
plt.plot(x,flux_new, label = '(input)')
plt.xlabel('wavelength (Angstrom)')
plt.ylabel('normalised flux (arb.)')
plt.legend()

#############################################################
#                                                           #
#                   Fourier Filter Plot                     #
#                                                           #
#############################################################
plt.figure()
x = np.linspace(0,500,500.0)
plt.scatter(x,abs(dft[:500]),s = 2)
plt.xlim(0,500)
plt.ylim(0,1e-13)
plt.ylabel('ampitude')
plt.xlabel('wavenumbers (k)')
plt.vlines(50,ymin=0,ymax=1e-13, linestyles='dashed', label='wavenumber cutoff')
plt.legend()
plt.savefig('fourierfreq.png')
plt.figure()

#########################################################
#                                                       #
#                   Correlation Plot                    #
#                                                       #
#########################################################

ax = plt.subplot(111)
l=[0,1024,2048]
x=["No Overlap","Full Overlap","No Overlap"]
ax.plot(correlation, label = 'Correlation')
ax.plot(a_n, label = 'antisymmetry function')
ax.set_xticks(l)
ax.set_xticklabels(x)
ax.set_ylim(-0.25,1.3)
#ax.vlines(1024,-1,2, linestyles = 'dashed')
ax.set_xlim(0,2047)
plt.legend()
plt.ylabel('normalisation correlation amplitude')
'''