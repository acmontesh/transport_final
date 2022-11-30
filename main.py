import numpy as np

from functions import *

# Define grid size
irows = 500
jcols = 100

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

# Initialize Temperature array
temp_array = initialize_temp_array(irows, jcols, form_temp_array)
lambda_sor = 0.5

max_error = 1000
error_old = 1001
i = 0
while (max_error>25) and (i < 2000) and (max_error < error_old):
    i += 1
    error_old = max_error
    # Set boundaries
    temp_array_old = np.copy(temp_array) # for error computing
    temp_array = set_temp_bc(temp_array, form_temp_array, r_array, vz_array, dr, dz, pipe_j, shoe_i, t_surf, q_top,
                             q_right, q0, q_left, k_mud, jcols, irows, alpha, lambda_sor)
    temp_array = comp_temp_ij(temp_array, vz_array, r_array, dr, dz, alpha, irows, lambda_sor)
    error_array = np.absolute(temp_array - temp_array_old)
    max_error = error_array.max()
    i_steps_new = i
print(temp_array)