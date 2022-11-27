##################################################################
#   SET OF FUNCTIONS TO SOLVE ENERGY EQ. IN FINAL PROJECT TP. ####
#############  VORTICITY TEAM ####################################
##################################################################
import numpy as np
def vz(r,inside=1):
    # vz is the z-velocity as a function of r in ft/s. 
    #r is the distance from the center of the drill pipe, in inches. 
    #inside is a boolean buffer that indicates whether the node is inside or outside of the DP. 
    import math
    if inside==1:
        return 0.32*(1-(r/2.52)**2.)
    else:
        return 175.7*(1-(r/3.25)**2.-1.811*math.log(3.25/r)) 



def thDiffusivity(k, rho, cpHat):
    #thDiffusivity is the thermal diffusivity in ft^2/s
    #k is the thermal conductivity in BTU/h-ft-degF
    #rho is the mud density in ppg (Lbm/gal)
    #cpHat is the mud heat capacity per mass unit in BTU/lbm-degF
    return 3.7E-5*k/(rho*cpHat)

def lambda1(vz, dr, alpha):
    #lambda1 is a linear term that multiplies part of the discretized version of PDE. 
    #lambda1 is dimensionless. 
    #vz is fluid velocity in ft/s
    #dr is the cell length in the grid, in ft. 
    #alpha is the mud thermal diffusivity in ft^2/s.
    return vz*dr/(4*alpha)

def lambda2(r, dr):
    #lambda2 is a linear term that multiplies part of the discretized version of PDE. 
    #lambda2 is dimensionless. 
    #r is the distance from the center of the drill pipe, in inches. 
    #dr is the cell length in the grid, in ft. 
    return dr*3/r

def comp_temp_ij(vz, r, dr, alpha, t_iplus1, t_iminus1, t_jplus1, t_jminus1):
    '''
    Returns Temperature at current index ij.
    :param vz: fluid velocity in ft/s
    :param r: distance from the center of the drill pipe, in inches
    :param dr: cell length in the grid, in ft
    :param alpha: mud thermal diffusivity in ft^2/s
    :param t_iplus1: Temperature in i+1, j
    :param t_iminus1: Temperature in i-1, j
    :param t_jplus1: Temperature in i, j+1
    :param t_jminus1: Temperature in i, j-1
    :return: t_ij (Temperature at current index ij)
    '''

    # Compute current lambdas
    lambda1_ij = lambda1(vz, dr, alpha)
    lambda2_ij = lambda2(r, dr)
    # Compute current Temperature
    t_ij = lambda1_ij(t_iplus1-t_iminus1) + 0.5*(t_iplus1+t_iminus1) - lambda2_ij(t_jplus1 - t_jminus1)
    return t_ij

def r_array(rmax, jcols):
    '''
    Return array in r direction
    :param rmax: Maximum radius
    :param jcols: Number of columns in r direction
    :return: r_array
    '''
#     rstep = np.linspace(0, 10, 5)
#     rstep_array = np.tile(rstep, (4, 1))
    return np.linspace(0, rmax, jcols)

def z_array(zmax, irows):
    '''
    Return array in z direction
    :param zmax: Maximum depth (max. z)
    :param irows: Number of rows in z direction
    :return: z_array
    '''
    return np.linspace(0, zmax, irows)

def initialize_temp_array(irows, jcols):
    '''
    Returns array of ones for initial temperature array
    :param irows: Number of rows in array
    :param jcols: Number of columns in array
    :return: temp_array
    '''
    return np.ones((irows, jcols))

