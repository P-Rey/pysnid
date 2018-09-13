# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 09:45:37 2018

@author: Peter
"""

import numpy as np
from scipy import signal, fftpack

def correlate_spec(log_input_spectrum, log_tempalte_spectrum):
    correlation = signal.correlate(log_input_spectrum, log_tempalte_spectrum, mode = 'full')
    correlation = correlation / (np.std(fftpack.fft(log_input_spectrum)) * np.std(fftpack.fft(log_input_spectrum)) )

    a_n = np.zeros(len(correlation))

    for j in range(1, int(len(correlation))):
        a_n[j] = correlation[-j] - correlation[j-1]

    rmsA = np.std(a_n[:int(len(a_n)/2)])

    h = max( correlation )
    r = h / ( np.sqrt(2) * rmsA )
    lap = np.log(log_input_spectrum[-1]/log_input_spectrum[0])
    rlap =  r * lap
    return correlation, rlap

'''
plt.plot(autocorr)
plt.plot(correlation)
print('rmsA', rmsA)
print('argmax',np.argmax(correlation))
shift = abs((2047 /2 )- np.argmax(correlation))
print('shift', shift)


#index_peak = np.argwhere(h)
middleIndex = (len(correlation) - 1)/2
index_peak = np.argmax(correlation) 
wave_bin = (x[-1] - x[0]) / num
middleIndex = (len(correlation) - 1)/2
bin_shift = middleIndex - index_peak
delta = wave_bin * bin_shift
#dwlog = np.log(x[-1]/x[0]) / len(x)

z = (np.exp( delta*dwlog )-1)
#z = (lamba_obs -lamba_act) / lambda_act 
'''


#def calc_redshift_from_crosscorr(crossCorr, nw, dwlog):
#    deltaPeak = np.argmax(crossCorr)
#    print("deltaPeak",deltaPeak)
    # z = np.exp(deltaPeak * dwlog) - 1 #equation 13 of Blondin)
#    zAxisIndex = np.concatenate((np.arange(-nw / 2, 0), np.arange(0, nw / 2)))
#    print("zAxisIndex ",zAxisIndex )
#    plt.figure()
#    plt.plot(zAxisIndex)
#    print("nw",nw/2)
#
#    if deltaPeak < nw / 2:
#        z = (np.exp(abs(zAxisIndex) * dwlog) - 1)[deltaPeak]
#    else:
#        z = -(np.exp(abs(zAxisIndex) * dwlog) - 1)[deltaPeak]
#
#    return z
#dwlog = np.log(wavep01[-1]/wavep01[0]) / len(wavep01)
#redshift = calc_redshift_from_crosscorr(correlation, len(log_wave), dwlog)
