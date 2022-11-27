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
q0 = 926 # BTU/h ft2
# Create input arrays
r_array = r_array(rmax, jcols)
z_array = z_array(zmax, irows)
pipe_j = search_index(r_array, r_dp) # j index for r = r_dp
shoe_i = search_index(z_array, z_shoe) # i index for z = z_shoe
form_temp_array = gen_form_temp_array(temp_grad, t_surf, z_array)
# Initialize Temperature array
temp_array = initialize_temp_array(irows, jcols)

