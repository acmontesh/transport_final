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
q0 = 926 # BTU/h ft2
# Create arrays
r_array = r_array(rmax, jcols)
z_array = z_array(zmax, irows)
pipe_j = r_array.searchsorted(r_dp, 'right') - 1 # j index for r = r_dp
# Initialize Temperature array
temp_array = initialize_temp_array(irows, jcols)

