from ReadFile import Read
#This imports the function Read from Readfile
from ReadFile import data
#This imports the variable data
import numpy as np
import astropy.units as u
'''The following function finds the value of x, index and xnew'''
def ParticleInfo(filename, particle_type, partical_number):
    x = data['x']
    #is the x value
    index = np.where(data['x']>2)
    #indexes the values to find particle
    xnew = data['x'][index]
    #finds the new x value
