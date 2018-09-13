# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 17:32:38 2018

@author: Peter
"""
import numpy as np
import argparse

def readspec(filename):
    # Read in spectrum
    spectrum = []
    with open(filename) as file:
        for line in file:
            x = line.split()
            spectrum.append([float(x[i]) for i in range(len(x))])
    # Convert to numpy array
    spectrum = np.array(spectrum)
    return spectrum

def parsereadspec():
    # Read in spectrum
    parser = argparse.ArgumentParser()
    parser.add_argument("input_spec")
    args = parser.parse_args()
    spectrum = []
    with open(args.input_spec) as file:
        for line in file:
            x = line.split()
            spectrum.append([float(x[i]) for i in range(len(x))])
    # Convert to numpy array
    spectrum = np.array(spectrum)
    return spectrum