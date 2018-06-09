# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:52:50 2018

@author: Peter
"""

def red_boi_loop(temp_wave, shift):
    
    #############################################################
    #                                                           #
    #   returns a new array linearly shifted based on inputs.   #
    #   0.005 is arbitrary and can be changedb                  #
    #                                                           #
    #############################################################
    
    temp_wave_shift = [element + shift*0.005 for element in temp_wave]

    return temp_wave_shift

'''
def get_z_value(data_wave, filteredflux_data, temp_wave, filteredflux_temp):

#    look into the integer indexing of py it might be fucking with the h value 
    
    r_list=[]
    N = 30          #THis means that we will look for z values up to z=3 over 0.1 increments 
    for i in range(0,N):
        print(i)
        temp_wave = red_boi_loop(temp_wave,N)
        r = Correlate(data_wave, filteredflux_data, temp_wave, filteredflux_temp)
        np.append(r_list,r)
    return rs
'''