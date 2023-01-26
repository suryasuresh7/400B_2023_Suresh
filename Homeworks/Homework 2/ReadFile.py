# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 20:45:37 2022

@author: surya
"""
import numpy as np
import astropy.units as u

def Read(filename):
    file = open(filename,'r' )
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    line2 = file.readline()
    label1, value1 = line2.split()
    total = float(value1)
    file.close()
    global data
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, total, data
Read("MW_000.txt")