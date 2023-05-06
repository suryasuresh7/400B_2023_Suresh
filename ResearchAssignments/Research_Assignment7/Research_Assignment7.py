
'''
This code, calculates the Merger of the MW and M31 and finds the stellar density
profile of the merger remnant. The main function, "plot_surface_density," 
loads particle data for the merger remnant of the Milky Way and Andromeda galaxies. 
It computes the surface density profiles for both galaxies and plots them along 
with the simulated bulge. Additionally, the code includes a function, 
"density_contour," to create a density contour plot. Finally, it produces a 
plot of the density distribution for the merger remnant using a 2D histogram 
and contour lines.'''
#importing modules
#importing modules
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#Code done in homeworks and class
from ReadFile import Read
from MassProfile_Soln import MassProfile
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass
#Function to calculate Sersic Profile
'''Taken from Lab 6'''
def sersicE(r, re, n, mtot):
    """ Function that computes the Sersic Profile for an Elliptical 
    System, assuming M/L ~ 1
    PARMETERS
    ---------
        r: `float`
            Distance from the center of the galaxy (kpc)
        re: `float`
            The Effective radius (2D radius that contains 
            half the light) (kpc)
        n:  `float`
            sersic index
        mtot: `float`
            the total stellar mass (Msun)

    RETURNS
    -------
        I: `array of floats`
            the surface brightness profile of the elliptical in Lsun/kpc^2

    """

    # We are assuming M/L = 1, so the total luminosity is:
    lum = mtot
    
    # the effective surface brightness is
    Ie = lum/7.2/np.pi/re**2
    
    # Break down the equation 
    a = (r/re)**(1.0/n)
    b = -7.67*(a-1)
    
    # The surface brightness
    #I = Ie*np.exp(-7.67*((r/re)**(1.0/n)-1.0))
    I = Ie*np.exp(b)
    
    return I
'''Original function, Lab 6 code modified and merger remnant simulation.'''
def plot_surface_density(filename):
    '''
    Function to plot the surface density profile of the MW+M31 merger remnant bulge.
    
    PARAMETERS
    ----------
    filename: `str`
        The filename used to load the particle data.
    
    RETURNS
    -------
    Plot of Stellar Density Profile
    '''
    # Load particle data for M31 and MW
    COM_M31 = CenterOfMass(f"M31_00{filename}.txt", 3)  # Center of Mass of M31
    COM_MW = CenterOfMass(f"MW_00{filename}.txt", 3)  # Center of Mass of MW
    
    # Compute the Center of Mass positions
    COM_M31_P = COM_M31.COM_P(0.1, 2.0)  # Center of Mass position of M31
    COM_MW_P = COM_MW.COM_P(0.1, 2.0)  # Center of Mass position of MW

    # Shift the coordinates to the center of mass positions
    x_M31 = COM_M31.x - COM_M31_P[0].value
    y_M31 = COM_M31.y - COM_M31_P[1].value
    z_M31 = COM_M31.z - COM_M31_P[2].value
    m_M31 = COM_M31.m  # Mass of M31

    x_MW = COM_MW.x - COM_MW_P[0].value
    y_MW = COM_MW.y - COM_MW_P[1].value
    z_MW = COM_MW.z - COM_MW_P[2].value
    m_MW = COM_MW.m  # Mass of MW

    # Combine the mass and position arrays of M31 and MW
    M_MW_M31 = np.append(m_MW, m_M31)
    x_MW_M31 = np.append(x_MW, x_M31)
    y_MW_M31 = np.append(y_MW, y_M31)
    z_MW_M31 = np.append(z_MW, z_M31)

    # Compute the cylindrical coordinates and radii
    cyl_r_mag = np.sqrt(x_MW_M31**2 + y_MW_M31**2)
    cyl_theta = np.arctan2(y_MW_M31, x_MW_M31)
    radii = np.arange(0.1, 0.95 * cyl_r_mag.max(), 0.1)
    
    # Compute the surface density profile of M31
    cyl_r_mag_M31 = np.sqrt(x_M31**2 + y_M31**2)
    enc_mask_M31 = cyl_r_mag_M31[:, np.newaxis] < np.asarray(radii).flatten()
    m_enc_M31 = np.sum(m_M31[:, np.newaxis] * enc_mask_M31, axis=0)
    m_annuli_M31 = np.diff(m_enc_M31)
    Sigma_M31 = m_annuli_M31 / (np.pi * (radii[1:] ** 2 - radii[:-1] ** 2))
    
    # Compute the surface density profile of MW
    cyl_r_mag_MW = np.sqrt(x_MW**2 + y_MW**2)
    enc_mask_MW = cyl_r_mag_MW[:, np.newaxis] < np.asarray(radii).flatten()
    m_enc_MW = np.sum(m_MW[:, np.newaxis] * enc_mask_MW, axis=0)
    m_annuli_MW = np.diff(m_enc_MW)
    Sigma_MW = m_annuli_MW / (np.pi * (radii[1:] ** 2 - radii[:-1] ** 2))

    # create the mask to select particles for each radius
    enc_mask = cyl_r_mag[:, np.newaxis] < np.asarray(radii).flatten()

    # calculate the enclosed masses
    m_enc = np.sum(M_MW_M31[:, np.newaxis] * enc_mask, axis=0)

    # use the difference between nearby elements to get mass in each annulus
    m_annuli = np.diff(m_enc)
    Sigma = m_annuli / (np.pi * (radii[1:] ** 2 - radii[:-1] ** 2))

    r_annuli = np.sqrt(radii[1:] * radii[:-1])

    # Plot surface density profile
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.loglog(r_annuli, Sigma, lw=2, alpha=0.8, label='Simulated Bulge')
    # Plot surface density profile of M31
    ax.loglog(r_annuli, Sigma_M31, lw=2, alpha=0.8, label='M31')

    # Plot surface density profile of MW
    ax.loglog(r_annuli, Sigma_MW, lw=2, alpha=0.8, label='MW')
    ax.set(xlabel=r"$r$ [kpc]", ylabel=r"$\Sigma$(disk par) [$10^{10} M_\odot$ / kpc$^2$]", 
           title=f"MW+M31 merger remnant bulge ({filename})")

    plt.xlim(0.1, 50)
    plt.ylim(1e-6, 0.2e1)

    ax.legend(loc='best')
    fig.tight_layout()

    plt.show()
    
plot_surface_density("0")

'''The below code has been taken from lab 7'''
#Now let's find the shape of the merger remnant.
# Code for plotting contours
#from https://gist.github.com/adrn/3993992

import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

# Info about **kwargs, *args 
#https://book.pythontips.com/en/latest/args_and_kwargs.html

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        i.e. unknown number of keywords 
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
             colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))  
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T  # transpose of distribution fxn
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    #brentq is root finding method
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    
    # You might need to add a few levels
    # I added a few between 1 and 2 sigma to better highlight the spiral arm
    one_sigma1 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.80))
    one_sigma2 = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.90))

    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, one_sigma1, one_sigma2, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.68', '0.8','0.9','0.95', '0.99'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour
'''The below is original'''
COM_M31 = CenterOfMass(f"M31_000.txt", 2)
COM_MW = CenterOfMass(f"MW_000.txt", 2)
COM_M31_P = COM_M31.COM_P(0.1, 2.0)
COM_MW_P = COM_MW.COM_P(0.1, 2.0)
x_M31 = COM_M31.x - COM_M31_P[0].value
y_M31 = COM_M31.y - COM_M31_P[1].value
z_M31 = COM_M31.z - COM_M31_P[2].value
m_M31 = COM_M31.m
x_MW = COM_MW.x - COM_MW_P[0].value
y_MW = COM_MW.y - COM_MW_P[1].value
z_MW = COM_MW.z - COM_MW_P[2].value
m_MW = COM_MW.m

M_MW_M31 = np.append(m_MW, m_M31)
x_MW_M31 = np.append(x_MW, x_M31)
y_MW_M31 = np.append(y_MW, y_M31)
z_MW_M31 = np.append(z_MW, z_M31)

COM_M31_V = COM_M31.COM_V(COM_M31_P[0],COM_M31_P[1],COM_M31_P[2])
COM_MW_V = COM_MW.COM_V(COM_MW_P[0],COM_MW_P[1],COM_MW_P[2])
vx_M31 = COM_M31.vx - COM_M31_V[0].value
vy_M31 = COM_M31.vy - COM_M31_V[1].value
vz_M31 = COM_M31.vz - COM_M31_P[2].value
vx_MW = COM_MW.vx - COM_MW_V[0].value
vy_MW = COM_MW.vy - COM_MW_V[1].value
vz_MW = COM_MW.vz - COM_MW_V[2].value

vx_MW_M31 = np.append(vx_MW, vx_M31)
vy_MW_M31 = np.append(vy_MW, vy_M31)
vz_MW_M31 = np.append(vz_MW, vz_M31)

'''The below code has been taken from Lab 7'''
# 1) Make plots 

# MW+M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# ADD HERE
# plot the particle density for MW+M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap= )
# cmap options: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. magma, viridis
# can modify bin number to make the plot smoother
plt.hist2d(x_MW_M31,y_MW_M31, bins= 1500, norm = LogNorm(), cmap = 'magma')

plt.colorbar()

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(x_MW_M31,y_MW_M31, 80, 80, ax=ax, colors =['red','orange','yellow','green'])


# Add axis labels
plt.xlabel(' ', fontsize=22)
plt.ylabel(' ', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size



# Save to a file
plt.savefig('Lab7_M31Disk.png')

#Rotating galaxy
def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel
# Vectors for r and v 
r = np.array([x_MW_M31,y_MW_M31,z_MW_M31]).T # transposed 
v = np.array([vx_MW_M31,vy_MW_M31,vz_MW_M31]).T
rn, vn = RotateFrame(r,v)
# MW+M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the particle density for MW+M31 
# ADD HERE
plt.hist2d(rn[:,0], rn[:,1], bins = 1500, norm = LogNorm(), cmap = 'viridis')
plt.colorbar()

# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# ADD HERE

# Add axis labels
plt.xlabel('  ', fontsize=22)
plt.ylabel('  ', fontsize=22)

#set axis limits
plt.ylim(-40,40)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size

# Save to a file 
plt.savefig('Lab7_FaceOn_Density.png')

# Rotated MW+M31 Disk - EDGE ON

# MW+M31 Disk Density 
fig, ax= plt.subplots(figsize=(10, 10))

# plot the w , 2D histogram
# ADD HERE
plt.hist2d(rn[:,0], rn[:,2], bins = 1500, norm = LogNorm(), cmap = 'viridis')

plt.colorbar()

# Add axis labels
plt.xlabel(' ', fontsize=22)
plt.ylabel(' ', fontsize=22)

#set axis limits
plt.ylim(-10,10)
plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size

# Save to a file
plt.savefig('Lab7_EdgeOn_Density.png')