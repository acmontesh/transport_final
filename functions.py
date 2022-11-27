##################################################################
#   SET OF FUNCTIONS TO SOLVE ENERGY EQ. IN FINAL PROJECT TP. ####
#############  VORTICITY TEAM ####################################
##################################################################

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