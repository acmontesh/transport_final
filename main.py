from functions import *

# Define grid size
irows = 10
jcols = 10

# Define well dimensions
rmax = 12.25/2 # inches
r_dp = 6/2 # inches
zmax = 2000 # ft
z_shoe = 1000 # ft

# Boundary Conditions
t_surf = 25 # degC
temp_grad = 3/3.28 # degC/100 ft
q0 = 926 # BTU/h ft2 - From pipe friction
q_top = 0
q_right = 0 # Casing isolation
q_left = 0

# Create input arrays
r_array, dr = r_array(rmax, jcols)
z_array, dz = z_array(zmax, irows)
pipe_j = search_index(r_array, r_dp) # j index for r = r_dp
shoe_i = search_index(z_array, z_shoe) # i index for z = z_shoe
form_temp_array = gen_form_temp_array(temp_grad, t_surf, z_array)
vz_array = comp_vz_array(r_array, pipe_j, irows, jcols)
# Initialize Temperature array
temp_array = initialize_temp_array(irows, jcols)

# Set boundaries
temp_array = set_temp_bc(temp_array, form_temp_array, dr, dz, pipe_j, shoe_i, t_surf, q_top, q_right, q0, q_left)
# # Upper Boundary
# temp_array[0, 0:pipe_j+1] = t_surf
# temp_array[0, pipe_j+1:] = temp_array[2, pipe_j+1:] - 2*dr*q_top
# # Right Boundary
# temp_array[:shoe_i+1, -1] = temp_array[:shoe_i+1, -3] - 2*dr*q_right
# temp_array[shoe_i+1:, -1] = form_temp_array[shoe_i+1:]
# # Bottom Boundary
# temp_array[-1, :] = temp_array[-3, :] - 2*dz*q0
# # Left Boundary
# temp_array[:, 0] = temp_array[:, 2] - 2*dr*q_left
print(temp_array)
