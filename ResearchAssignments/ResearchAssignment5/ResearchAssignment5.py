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
from MassProfile_Soln import MassProfile
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass
#Function to calculate Sersic Profile
class TotalMass:
    #May need other inputs not sure yet
    def __init__(self, snap_number_M31, snap_number_MW):
        # Create instances for MassProfile class for M31 and MW
        self.Mass_M31 = MassProfile("M31", snap_number_M31)
        self.Mass_MW = MassProfile("MW", snap_number_MW)

        # Read data for M31
        self.time, self.total, self.data = Read("M31_000.txt")
        self.M_31 = self.data['m'] * u.Msun
        self.x_31 = self.data['x'] * u.kpc
        self.y_31 = self.data['y'] * u.kpc
        self.z_31 = self.data['z'] * u.kpc

        # Read data for MW
        self.time_MW, self.total_MW, self.data_MW = Read("MW_000.txt")
        self.M_MW = self.data_MW['m'] * u.Msun
        self.x_MW = self.data_MW['x'] * u.kpc
        self.y_MW = self.data_MW['y'] * u.kpc
        self.z_MW = self.data_MW['z'] * u.kpc

        # Combine particle data of M31 and MW to simulate merger
        self.M_MW_M31 = np.append(self.M_31, self.M_MW)
        self.x_MW_M31 = np.append(self.x_MW, self.x_31)
        self.y_MW_M31 = np.append(self.y_MW, self.y_31)
        self.z_MW_M31 = np.append(self.z_MW, self.z_31)

    def Massenclosed(self, ptype, r_0, r_f,n):
        com = CenterOfMass("M31_000.txt",2)
        com_pos = com.COM_P(0.1, 2.0)
        index = np.where(self.data['type'] == ptype)
        #finding CenterOfMass for remnant.
        xG = self.x_MW_M31[index] - com_pos[0]
        yG = self.y_MW_M31[index] - com_pos[1]
        zG = self.z_MW_M31[index] - com_pos[2]
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
        mG = self.M_MW_M31[index]
        mG = mG.value
        #Let radii increase by a small amount every iteration (delta r) 
        #from 0.05 kpc to 60 kpc. Thus giving a small change instead of 
        #taking the whole radius.
        radii = np.linspace(r_0,r_f,n)
        m_enc = np.zeros(np.size(radii))
        density_r = np.zeros(np.size(radii -1))
        for i in range(len(radii))[:-1]:
            indexR = np.where((rG >= radii[i]*u.kpc) & (rG <= radii[i+1]*u.kpc))
            m_enc[i] = np.sum(mG[indexR])
            shell_area = 4 * np.pi * radii[i] ** 2  # Surface area of the shell
            density = m_enc[i] / shell_area  # Density in the shell
            density_r[i+1] = density
            i+=1
        return (density_r * (1e10 * u.Msun/u.pc**3)) 
    def sersicE(r, re, n, mtot):
        #Total luminosity
        lum = mtot
        #Effective surface brightness is
        Ie = lum/7.2/np.pi/re**2
        a = (r/re)**(1.0/n)
        b = -7.67*(a-1)
        # The surface brightness
        I = Ie*np.exp(-7.67*((r/re)**(1.0/n)-1.0))
        I = Ie*np.exp(b)
        return I
    

totalmass = TotalMass(0,0)
r_0 = 0 #initial radius
r_f = 20 #final radius
n = 100 #number of radii
# Call the Massenclosed method with desired particle type
stellar_density = totalmass.Massenclosed(2, r_0, r_f, n)
# Plot the density profile
radii = np.linspace(r_0, r_f, n)
plt.plot(radii, stellar_density * 1e10, label='Stellar Density')
plt.xlabel('Radius (kpc)')
plt.ylabel('Density (Msun/pc^3)')
plt.title('Surface Density Profile (Computed in Shells)')
plt.legend()
plt.show()