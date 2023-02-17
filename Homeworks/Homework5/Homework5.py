# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:26:13 2023

@author: surya
"""

from ReadFile import Read
from CenterOfMass import CenterOfMass
import astropy.units as u
import astropy.constants as const
import numpy as np
import matplotlib.pyplot as plt
class MassProfile:
    '''This class contains all the function required to calculate the
    Mass enclosed and Circular velocity'''
    def __init__(self,galaxy,snap):
        '''Initializes the class, takes self, galaxy which is the name of 
        the galaxy, snap is the number after the Galaxy name for example"0".
        '''
        # add a string of the filenumber to the value “000”
        ilbl = "000" + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename="%s "%(galaxy) + ilbl + ".txt"
        self.time, self.total, self.data = Read(self.filename)
        self.m = self.data['m'] #gets the mass
        self.x = self.data['x']*u.kpc #gets the x pos
        self.y = self.data['y']*u.kpc #gets the y pos
        self.z = self.data['z']*u.kpc #gets the z pos
        self.vx = self.data['vx']*u.km/u.s #gets the x velocity
        self.vy = self.data['vy']*u.km/u.s #gets the y velocity
        self.vz = self.data['vz']*u.km/u.s #gets the z velocity
        self.gname = galaxy #the name of galaxy
    def MassEnclosed(self,ptype, radii):
        '''This function takes the inputs:
            ptype = the particle type (0 = bulge, 1 = disk, 2 = halo)
            radii = the radii is the list of the radius
            returns the mass enclosed in units of MSun
        '''
        MW_COM = CenterOfMass(self.filename,ptype) 
        MW_COM1 = MW_COM.COM_P(0.1)
        mass = np.zeros(len(radii))
        r = np.sqrt((self.x - MW_COM1[0])**2 + (self.y - MW_COM1\
        [1])**2 + (self.z - MW_COM[2])**2)
        for i in range(len(radii)): #finds the mass enclosed
            idx = np.where(r < radii[i])
            m = np.sum(self.m[idx])
            mass[i] = m
        return mass*u.Msun 
    def MassEnclosedTotal(self,radii):
        '''This function takes the inputs:
            radii = the radii is the list of the radius
            returns the total mass enclosed, the total mass enclosed in
            bulge, the total mass enclosed in disk, and halo.
        '''
        m_bulge = self.MassEnclosed(1, radii) #the mass of bulge
        m_disk = self.MassEnclosed(2, radii) #the mass of disk
        m_halo = self.MassEnclosed(3, radii) #the mass of halo
        if self.gname == "M33": 
            m_total = m_disk + m_halo 
        else:
            m_total = np.sum(m_bulge, m_disk, m_halo)
        return m_total, m_bulge, m_disk, m_halo
    def HernquistMass(self, radii,a, M_halo):
        '''This function takes the inputs:
            radii = the radii is the list of the radius
            a = scale factor
            M_halo = mass of halo
            Finds the Hernquist mass radially'''
        mass_r = (M_halo*(radii**2))/(a+radii)**2
        return mass_r
    def CircularVelocity(self, ptype, radii):
        '''Finds the Circular velocity in km/s
            This function takes the inputs:
            ptype = the particle type (0 = bulge, 1 = disk, 2 = halo)
            radii = the radii is the list of the radius
            returns v_circ (circular velocity)
        '''
        mass_1 = self.MassEnclosed(ptype,radii)
        v_circ = np.zeros(len(radii))
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) #G constant 
        for i in range(len(radii)):
            m_1 = np.sum(mass_1[i]) #If not work use self.m[i]
            v = np.sqrt(G * m_1/radii[i]) #velocity formula
            v_circ[i] = v
        return v_circ*u.km/u.s
    def CircularVelocityTotal(self,radii):
        '''Finds the total circular velocity in km/s
        This function takes the inputs:
        radii = the radii is the list of the radius
        returns v_circ1 (total circular velocity)
        '''
        v_circ = np.zeros(len(radii))
        v_circ1 = np.zeros(len(radii))
        v_circ2 = np.zeros(len(radii))
        mass_bulge = self.MassEnclosedTotal(radii)[1] #mass of bulge
        mass_disk = self.MassEnclosedTotal(radii)[2] #mass of disk
        mass_halo = self.MassEnclosedTotal(radii)[3] #nass of halo
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) #G constant 
        for i in range(len(radii)):
            v_b = np.sqrt(G * mass_bulge/radii[i]) #velocity of bugle
            v_d = np.sqrt(G * mass_disk/radii[i]) #velocity of disk
            v_h = np.sqrt(G * mass_halo/radii[i]) #velocity of halo
            v_circ1 = [v_b,v_d,v_h] #circular velocity
        return v_circ1*u.km/u.s
    def HernquistVCirc(self, radii, a ,M_halo):
        '''Finds the total Hernquist circular velocity in km/s
        This function takes the inputs:
        radii = the radii is the list of the radius
        returns v_circ1 (total circular velocity)
        '''
        Hernquist_mass = self.HernquistMass(radii,a, M_halo) 
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        v_circ = np.sqrt(G * Hernquist_mass/radii)

#Plot Mass Profile
radii = np.arange(0,30,20)
mass = MassProfile("MW", 0).MassEnclosed(1,radii)
plt.plot(radii, mass, 'red')


            
            
        
            
            