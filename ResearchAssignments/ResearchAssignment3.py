# -*- coding: utf-8 -*-

'''
Name: Surya Suresh
'''
"""
Stellar density Profiles of merger remnants.

The topic of my research Project is, "MW+M31 Stellar Major Merger Remnant: 
Stellar disk particle distribution/morphology." The question I have chosen is:
"What is the final stellar density profile for the combined system ? Is it well 
fit by a sersic profile? Does it agree with predictions for elliptical 
galaxies?"
For this topic, I intend to make a plot that plots the Stellar density profile 
of the merger remnant as a function of radius. For this plot, lab 6 provides 
a lot of guidance, but for this research project the mass profile and stellar 
density profile needs to be computed in shells rather than the whole bulge. 
"""
#importing modules
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
#Code done in homeworks and class
from ReadFile import Read
from MassProfile import MassProfile
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass
#Function to calculate Sersic Profile
def Sersic(r, re, n, M):
    """ This Function finds the Sersic Profile for a given galaxy, in our case
    it should calculate the sersic profile of the MW + M31 Merger Remnant. It 
    takes inputs r (radius that we are computing the sersic profile for), re 
    (the effective radius), n (the Sersic index) and M (The total Mass)."""
    Luminosity = M
    Ie = Luminosity/7.2/np.pi/re**2
    I_tot = Ie*np.exp(-7.67((r/re)**(1/n))-1)
    return I_tot
    #The issue with this approach is that this assumes that the density 
    #p = M/V is constant, as approached by the MassProfile code. The code 
    #needed is to compute density in shells. Let us try to recreate this in 
    #Mass profile code. This code must be edited to accomodate this restriction.
class TotalMass:
    def __init__(self,'''need inputs'''):
        '''Needs variables'''
        Mass_M31 = MassProfile("M31", '''needs snap number''')
        Mass_MW = MassProfile("MW", '''needs snap number''')
        #Need to figure out ways to combine particle data of M31 and Milky Way
        #Assume for now we add them like arrays.
        self.time, self.total, self.data = Read("M31_000.txt")
        #Instance for M31
        self.M_31 = self.data['m']*u.Msun
        self.x_31 = self.data['x']*u.kpc
        self.y_31 = self.data['y']*u.kpc
        self.z_31 = self.data['z']*u.kpc
        self.time_MW, self.total_MW, self.data_MW = Read("MW_000.txt")
        #Instance for MW
        self.M_MW = self.data_MW['m']*u.Msun
        self.x_MW = self.data_MW['x']*u.kpc
        self.y_MW = self.data_MW['y']*u.kpc
        self.z_MW = self.data_NW['z']*u.kpc
        #Adding arrays to get MW and M31 merger
        self.M_MW_M31 = np.append(self.M_31,self.M_MW)
        self.x_MW_M31 = np.append(self.x_MW,self.x_31)
        self.y_MW_M31 = np.append(self.y_MW,self.y_31)
        self.z_MW_M31 = np.append(self.z_MW,self.z_31)
    def Massenclosed(self, ptype):
        com = CenterOfMass("M31_000.txt",2)
        com_pos = com.COM_P(0.1)
        index = np.where(self.data['type'] == ptype)
        #finding CenterOfMass for remnant.
        xG = self.x_MW_M31[index] - com_pos[0]
        yG = self.y_MW_M31[index] - com_pos[1]
        zG = self.z_MW_M31[index] - com_pos[2]
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
        mG = self.M_MW_M31[index]
        #Have to figure out how to compute mass enclosed in shells.
        m_enc = np.zeros(np.size(radii))
        #Let radii increase by 0.01kpc every iteration (delta r) from 0.05 kpc
        #to 60 kpc. Thus giving a small change instead of taking the whole 
        #radius.
        radii = np.linspace(0.05,60,200)
        for i in range(len(radii))[:-1]:
            indexR = np.where(rG >= radii[i]*u.kpc & rG <= radii[i+1])
            m_enc[i] = np.sum(mG[indexR])
            i+=1
        return m_enc*u.Msun*1e10
        
        
