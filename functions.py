##################################################################
#   SET OF FUNCTIONS TO SOLVE ENERGY EQ. IN FINAL PROJECT TP. ####
#############  VORTICITY TEAM ####################################
##################################################################

def vz(r,inside=1, holeR=3.25, pipeR=2.52):
    import numpy as np
    # vz is the z-velocity as a function of r in ft/s. 
    #r is the distance from the center of the drill pipe, in inches. 
    #inside is a boolean buffer that indicates whether the node is inside or outside of the DP.
    # Velocity function based on MW 10 ppg assumption / 10 gpm flowrate
    # R pipe 2.52 - R well 3.25
    # import math
    if inside==1:
        return 0.32*(1-(r/pipeR)**2.)
    else:
        kappa = pipeR/holeR
        return 175.7*(1-(r/holeR)**2.-((1-kappa**2.)/np.log(1/kappa))*np.log(holeR/r))

def comp_vz_array(r_array, pipe_j, irows, jcols):
    import numpy as np
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
    vz_array[:,pipe_j+1:] = -vz(r_array[:,pipe_j+1:], inside=0)
    return vz_array

def thDiffusivity(k, rho, cpHat):
    #thDiffusivity is the thermal diffusivity in ft^2/s
    #k is the thermal conductivity in BTU/h-ft-degF
    #rho is the mud density in ppg (Lbm/gal)
    #cpHat is the mud heat capacity per mass unit in BTU/lbm-degF
    return 3.7E-5*k/(rho*cpHat)

def lambda1(vz, dr, dz, alpha):
    #lambda1 is a linear term that multiplies part of the discretized version of PDE. 
    #lambda1 is dimensionless. 
    #vz is fluid velocity in ft/s
    #dr is the cell length in the grid, in ft. 
    #alpha is the mud thermal diffusivity in ft^2/s.
    return vz*(dr**2)/(4*alpha*dz)

def lambda2(r, dr):
    #lambda2 is a linear term that multiplies part of the discretized version of PDE. 
    #lambda2 is dimensionless. 
    #r is the distance from the center of the drill pipe, in inches. 
    #dr is the cell length in the grid, in ft. 
    return dr*3/r

def comp_temp_ij(temp_array, vz_array, r_array, dr, dz,alpha, irows, lambda_sor=1.5):
    import numpy as np
    '''
    Returns Temperature Array 'temp_array' after computing heat transfer
    :param vz: fluid velocity in ft/s
    :param r: distance from the center of the drill pipe, in inches
    :param dr: cell length in the grid, in ft
    :param alpha: mud thermal diffusivity in ft^2/s
    :return: temp_array
    '''

    # Compute current lambdas
    lambda1_ij = lambda1(vz_array, dr, dz, alpha)
    lambda1_ij = lambda1_ij[1:-1, 1:-1].copy()
    lambda2_ij = lambda2(r_array, dr)
    lambda2_ij = np.tile(lambda2_ij, (irows, 1))
    lambda2_ij = lambda2_ij[1:-1, 1:-1].copy()
    # Compute current Temperature
    t_iplus1 = temp_array[2:, 1:-1]
    t_iminus1 = temp_array[:-2, 1:-1]
    t_jplus1 = temp_array[1:-1, 2:]
    t_jminus1 = temp_array[1:-1, :-2]

    temp_array_new = lambda2_ij*(t_iplus1 - t_iminus1) + 0.5 * (t_iplus1 + t_iminus1) - lambda1_ij*(
            t_jplus1 - t_jminus1)
    # SOR
    temp_array_sor = lambda_sor*temp_array_new + (1-lambda_sor)*temp_array[1:-1, 1:-1]
    temp_array[1:-1, 1:-1] = temp_array_sor
    return temp_array

def r_array(rmax, jcols):
    import numpy as np
    '''
    Return array in r direction, and dr
    :param rmax: Maximum radius
    :param jcols: Number of columns in r direction
    :return: r_array, dr
    '''

    r_array_i = np.linspace(0.0001, rmax, jcols)
    dr = rmax/(jcols - 1)
    return r_array_i, dr

def z_array(zmax, irows):
    import numpy as np
    '''
    Return array in z direction, and dz
    :param zmax: Maximum depth (max. z)
    :param irows: Number of rows in z direction
    :return: z_array, dz
    '''
    z_array_i = np.linspace(0, zmax, irows)
    dz = zmax/(irows - 1)
    return z_array_i, dz


def initialize_temp_array(irows, jcols, form_temp_array):
    import numpy as np
    '''
    Returns array of ones for initial temperature array
    :param irows: Number of rows in array
    :param jcols: Number of columns in array
    :return: temp_array
    '''
    temp_array = np.tile(form_temp_array.reshape(-1,1), (1, jcols))
    return temp_array

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


def set_temp_bc(temp_array, form_temp_array, dr, dz, pipe_j, shoe_i, t_surf, q_top, q_right, q0, q_left, k_mud):
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
    temp_array[0, pipe_j + 1:] = temp_array[2, pipe_j + 1:] - 2 * dz * q_top
    # Right Boundary
    temp_array[:shoe_i + 1, -1] = temp_array[:shoe_i + 1, -3] - 2 * dr * q_right
    temp_array[shoe_i + 1:, -1] = form_temp_array[shoe_i + 1:]
    # Bottom Boundary
    temp_array[-1, :] = temp_array[-3, :] + q0/k_mud # Based on Fourier's law
    # Left Boundary
    temp_array[:, 0] = temp_array[:, 2] - 2 * dr * q_left
    return temp_array