##################################################################
#   SET OF FUNCTIONS TO SOLVE ENERGY EQ. IN FINAL PROJECT TP. ####
#############  VORTICITY TEAM ####################################
##################################################################
import numpy as np
def vz(r,inside=1):
    # vz is the z-velocity as a function of r in ft/s. 
    #r is the distance from the center of the drill pipe, in inches. 
    #inside is a boolean buffer that indicates whether the node is inside or outside of the DP. 
    # import math
    if inside==1:
        return 0.32*(1-(r/2.52)**2.)
    else:
        return 175.7*(1-(r/3.25)**2.-1.811*np.log(3.25/r))

def comp_vz_array(r_array, pipe_j, irows, jcols):
    '''
    Return vz array
    :param r_array: Dimensions array in r direction
    :param pipe_j: Index for pipe wall
    :param irows: i-rows for grid dimension
    :param jcols: j-cols for grid dimension
    :return: vz_array
    '''

    r_array = np.tile(r_array, (irows, 1))
    vz_array = np.ones((irows, jcols))
    vz_array[:,:pipe_j+1] = vz(r_array[:,:pipe_j+1], inside=1)
    vz_array[:,pipe_j+1:] = vz(r_array[:,pipe_j+1:], inside=0)
    return vz_array

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

def comp_temp_ij(temp_array, vz, r, dr, alpha):
    '''
    Returns Temperature Array 'temp_array' after computing heat transfer
    :param vz: fluid velocity in ft/s
    :param r: distance from the center of the drill pipe, in inches
    :param dr: cell length in the grid, in ft
    :param alpha: mud thermal diffusivity in ft^2/s
    :return: temp_array
    '''

    # Compute current lambdas
    lambda1_ij = lambda1(vz, dr, alpha)
    lambda2_ij = lambda2(r, dr)
    # Compute current Temperature
    #t_ij = lambda1_ij(t_iplus1-t_iminus1) + 0.5*(t_iplus1+t_iminus1) - lambda2_ij(t_jplus1 - t_jminus1)
    t_iplus1 = temp_array[2:, 1:-1]
    t_iminus1 = temp_array[:-2, 1:-1]
    t_jplus1 = temp_array[1:-1, 2:]
    t_jminus1 = temp_array[1:-1, :-2]
    temp_array[1:-1, 1:-1] = lambda1_ij(t_iplus1-t_iminus1) + 0.5*(t_iplus1+t_iminus1) - lambda2_ij(t_jplus1 - t_jminus1)
    return temp_array

def r_array(rmax, jcols):
    '''
    Return array in r direction, and dr
    :param rmax: Maximum radius
    :param jcols: Number of columns in r direction
    :return: r_array, dr
    '''
#     rstep = np.linspace(0, 10, 5)
#     rstep_array = np.tile(rstep, (4, 1))
    r_array_i = np.linspace(0, rmax, jcols)
    dr = rmax/(jcols - 1)
    return r_array_i, dr

def z_array(zmax, irows):
    '''
    Return array in z direction, and dz
    :param zmax: Maximum depth (max. z)
    :param irows: Number of rows in z direction
    :return: z_array, dz
    '''
    z_array_i = np.linspace(0, zmax, irows)
    dz = zmax/(irows - 1)
    return z_array_i, dz


def initialize_temp_array(irows, jcols):
    '''
    Returns array of ones for initial temperature array
    :param irows: Number of rows in array
    :param jcols: Number of columns in array
    :return: temp_array
    '''
    return np.ones((irows, jcols))

def search_index(arr, value):
    '''
    Return index number for value in array arr, where values are lower than value
    :param arr: Search array
    :param value: Required value
    :return: idx
    '''
    idx = arr.searchsorted(value, 'right') - 1
    return idx

def gen_form_temp_array(temp_grad, t_surf, z_array):
    '''
    Return formation temperature array based on temperature gradient 'temp_grad'
    and surface temperature 't_surf'
    :param temp_grad: Temperature Gradient in degC/ft
    :param t_surf: Surface Temperature in degC
    :param z_array: z_array (depth array)
    :return:formation temperature array
    '''

    form_temp_array = z_array*temp_grad/100 + t_surf
    return form_temp_array


def set_temp_bc(temp_array, form_temp_array, dr, dz, pipe_j, shoe_i, t_surf, q_top, q_right, q0, q_left):
    '''
    Return temp_array with boundary conditions
    :param temp_array: Temperature Array
    :param form_temp_array: Formation Temperature Array
    :param dr: Step in r direction
    :param dz: Step in z direction
    :param pipe_j: Pipe wall j index (column)
    :param shoe_i: Casing Shoe i index (row)
    :param t_surf: Mud Surface Temperature
    :param q_top: Heat transfer coefficient for annular top
    :param q_right: Heat transfer coefficient for r = 0 (left boundary)
    :param q0: Heat transfer coefficient due to pipe friction
    :param q_left: Heat transfer coefficient for Casing (right boundary)
    :return:
    '''
    # Set boundaries
    # Upper Boundary
    temp_array[0, 0:pipe_j + 1] = t_surf
    temp_array[0, pipe_j + 1:] = temp_array[2, pipe_j + 1:] - 2 * dr * q_top
    # Right Boundary
    temp_array[:shoe_i + 1, -1] = temp_array[:shoe_i + 1, -3] - 2 * dr * q_right
    temp_array[shoe_i + 1:, -1] = form_temp_array[shoe_i + 1:]
    # Bottom Boundary
    temp_array[-1, :] = temp_array[-3, :] - 2 * dz * q0
    # Left Boundary
    temp_array[:, 0] = temp_array[:, 2] - 2 * dr * q_left
    return temp_array