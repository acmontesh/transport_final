import numpy as np
from functions import *

# Define grid size
irows = 100
jcols = 10

# Fluid Properties
k_mud = 1.2*0.5778 # 1.2 W/m.K to BTU/h.ft.degF - From Magdy Abdel Hafis
rho = 10 # ppg - Assumption in Velocity Formula
cpHat = 0.3824 # 1600 J/kg.K to BTU/lbm-degF (1 to 0.00024) - From Magdy Abdel Hafis
# Review unit conversion

# Define well dimensions
rmax = 3.25/12 #ft # 12.25/2 # inches
r_dp = 2.52/12 #ft # inches
zmax = 2000 # ft
z_shoe = 1000 # ft

# Boundary Conditions
t_surf = 77 # degC
temp_grad = 1.64592 # degF/100 ft (equivalent to 3 degC/100m)
q0 = 926 # BTU/h ft2 - From pipe friction
q_top = 0
q_right = 0 # Casing isolation
q_left = 0

# Create input arrays
r_array, dr = gen_r_array(rmax, jcols)
r_array_inch = r_array*12
z_array, dz = gen_z_array(zmax, irows)
pipe_j = search_index(r_array, r_dp) # j index for r = r_dp
shoe_i = search_index(z_array, z_shoe) # i index for z = z_shoe
form_temp_array = gen_form_temp_array(temp_grad, t_surf, z_array)
vz_array = comp_vz_array(r_array_inch, pipe_j, irows, jcols) # r input must be in inches to compute velocity
alpha = thDiffusivity(k_mud, rho, cpHat)

# Arrays for plot
r_array_plot = np.tile(r_array, irows)
z_array_plot = np.repeat(z_array, jcols)
# Initialize Temperature array
temp_array = initialize_temp_array(irows, jcols, form_temp_array)
lambda_sor = 1.5

# Set Boundary Conditions
temp_array = set_temp_bc(temp_array, form_temp_array, r_array, vz_array, dr, dz, pipe_j, shoe_i, t_surf, q_top,
                             q_right, q0, q_left, k_mud, jcols, irows, alpha, lambda_sor)

# Jacobian Grid
j_grid_size = irows*jcols
j_grid = np.zeros((j_grid_size, j_grid_size))
# Results Vector
b_vector = np.zeros((j_grid_size))

for i in np.arange(0, j_grid_size):
    for j in np.arange(0, j_grid_size):
        # irow and jcol for original arrays
        irow_number = j // jcols
        jcol_number = j % jcols
        vz = vz_array[irow_number, jcol_number]
        r = r_array[jcol_number]
        lambda1_ij = lambda1(vz, dr, dz, alpha)
        lambda2_ij = lambda2(r, dr)
        # Boundaries
        if (i < jcols) or (i % jcols == 0) or ((i+1) % jcols == 0) or (i > (jcols*(irows-1))):
            # Set Boundary Temp in b_vector vector
            if i == j:
                b_vector[i] = temp_array[irow_number, jcol_number]
            # Set values in Jacobian Matrix
            if j == t_iplus1_j(i, jcols):
                j_grid[i, j] = lambda2_ij + 0.5
            elif j == t_iminus1_j(i, jcols):
                j_grid[i, j] = -lambda2_ij + 0.5
            elif j == t_i_jplus1(i, jcols):
                j_grid[i, j] = -lambda1_ij
            elif j == t_i_jminus1(i, jcols):
                j_grid[i, j] = lambda1_ij
        # Inside Grid (not boundary)
        else:
            if j == t_ij(i):
                j_grid[i, j] = -1
            elif j == t_iplus1_j(i, jcols):
                j_grid[i, j] = lambda2_ij + 0.5
            elif j == t_iminus1_j(i, jcols):
                j_grid[i, j] = -lambda2_ij + 0.5
            elif j == t_i_jplus1(i, jcols):
                j_grid[i, j] = -lambda1_ij
            elif j == t_i_jminus1(i, jcols):
                j_grid[i, j] = lambda1_ij

temp_array = jacobi(j_grid, b_vector)

# # Compute temperature grid - point by point
# for i in range(1, irows-1, 1)[::-1]:
#     for j in range(1, jcols-1, 1)[::-1]:
#         r_ij = r_array[j]
#         vz_ij = vz_array[i, j]
#         temp_array = comp_temp_ij_for_loop(temp_array, vz_ij, r_ij, dr, dz, alpha, lambda_sor, i, j)
#
# # Array for plot
# temp_array_plot = temp_array.flatten()
# # Save Arryas for plot
# directory = r'C:\Users\SaANTIAGO\Google Drive Streaming\My Drive\17_UT_Austin_PGE\00_Classes\381M_Transport_Phenomena\04_Final_Project\01_Code\transport_final'
# temp_array_file = '\\'.join([directory, 'temp_array.npy'])
# r_array_file = '\\'.join([directory, 'r_array.npy'])
# z_array_file = '\\'.join([directory, 'z_array.npy'])
# vz_array_file = '\\'.join([directory, 'vz_array.npy'])
# np.save(temp_array_file, temp_array)
# np.save(r_array_file, r_array_inch)
# np.save(z_array_file, z_array)
# np.save(vz_array_file, vz_array)

print(temp_array)