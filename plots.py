import numpy as np
import matplotlib.pyplot as plt

def plot_fig(log_temp = False, save_fig = True):
    # Define well dimensions
    rmax = 3.25 / 12  # ft # 12.25/2 # inches
    r_dp = 2.52 / 12  # ft # inches
    zmax = 2000  # ft
    z_shoe = 1000  # ft

    # Save Arryas for plot
    directory = r'C:\Users\SaANTIAGO\Google Drive Streaming\My Drive\17_UT_Austin_PGE\00_Classes\381M_Transport_Phenomena\04_Final_Project\01_Code\transport_final'
    temp_array_file = '\\'.join([directory, 'temp_array.npy'])
    r_array_file = '\\'.join([directory, 'r_array.npy'])
    z_array_file = '\\'.join([directory, 'z_array.npy'])
    temp_array = np.load(temp_array_file)
    r_array_inch = np.load(r_array_file)
    z_array = np.load(z_array_file)

    # Plot
    fig = plt.figure(figsize=(9,6))
    plt.title("Wellbore Temperature - Steady State - 1000 x 200 Grid - SOR 1.5")
    # sns.scatterplot(r_array_plot, z_array_plot, hue = temp_array_plot)
    if log_temp == True:
        temp_array = np.log(temp_array)
    plt.pcolormesh(r_array_inch, z_array, temp_array,
                   shading='gouraud',
                   cmap= "RdYlBu_r") #'plasma')

    # Plot Casing
    plt.plot([r_dp*12, rmax*12], [z_shoe, z_shoe], '--', color='gray')
    plt.plot([rmax*12, rmax*12], [0, z_shoe], '-', color='gray', lw=5)
    plt.text(rmax*12*.97, z_shoe*.8, 'Casing', rotation=90, color='gray')

    # Plot Drillpipe
    plt.plot([r_dp*12, r_dp*12], [0, zmax], '-', color='white', lw=3)
    plt.text(r_dp*12*.97, zmax/2, 'Drillpipe', rotation=90, color='white')

    # Colorbar
    cbar = plt.colorbar()
    cbar_label = 'Mud Temperature [degF]'
    if log_temp == True:
        cbar_label = 'Log Mud Temperature [degF]'
    cbar.ax.set_ylabel(cbar_label, rotation=90)

    # Axis
    plt.gca().invert_yaxis()
    plt.ylabel('Depth [ft]')
    plt.xlabel('Well Radius [in]')

    fig_path = '\\'.join([directory, 'temp_plot_1000 x 200 grid_SOR 1-5.png'])
    if save_fig == True:
        plt.savefig(fig_path)
    plt.show()
    return

def plot_vel(log_temp = False, save_fig = True):
    # Define well dimensions
    rmax = 3.25 / 12  # ft # 12.25/2 # inches
    r_dp = 2.52 / 12  # ft # inches
    zmax = 2000  # ft
    z_shoe = 1000  # ft

    # Save Arryas for plot
    directory = r'C:\Users\SaANTIAGO\Google Drive Streaming\My Drive\17_UT_Austin_PGE\00_Classes\381M_Transport_Phenomena\04_Final_Project\01_Code\transport_final'
    vz_array_file = '\\'.join([directory, 'vz_array.npy'])
    r_array_file = '\\'.join([directory, 'r_array.npy'])
    z_array_file = '\\'.join([directory, 'z_array.npy'])
    vz_array = np.load(vz_array_file)
    r_array_inch = np.load(r_array_file)
    z_array = np.load(z_array_file)

    # Plot
    fig = plt.figure(figsize=(9,6))
    plt.title("Fluid Velocity in z-direction ($V_z$) - Profile in r-direction - Steady State")
    # sns.scatterplot(r_array_plot, z_array_plot, hue = temp_array_plot)
    if log_temp == True:
        vz_array = np.log(vz_array)
    plt.pcolormesh(r_array_inch, z_array, vz_array,
                   shading='gouraud',
                   cmap= "RdYlBu_r") #'plasma')

    # Plot Casing
    plt.plot([r_dp*12, rmax*12], [z_shoe, z_shoe], '--', color='gray')
    plt.plot([rmax*12, rmax*12], [0, z_shoe], '-', color='gray', lw=5)
    plt.text(rmax*12*.97, z_shoe*.8, 'Casing', rotation=90, color='gray')

    # Plot Drillpipe
    plt.plot([r_dp*12, r_dp*12], [0, zmax], '-', color='white', lw=3)
    plt.text(r_dp*12*.97, zmax/2, 'Drillpipe', rotation=90, color='white')

    # Colorbar
    cbar = plt.colorbar()
    cbar_label = 'Fluid Velocity $V_z$ [ft/s]'
    if log_temp == True:
        cbar_label = 'Log Mud Temperature [degF]'
    cbar.ax.set_ylabel(cbar_label, rotation=90)

    # Axis
    plt.gca().invert_yaxis()
    plt.ylabel('Depth [ft]')
    plt.xlabel('Well Radius [in]')

    fig_path = '\\'.join([directory, 'velocity_plot.png'])
    if save_fig == True:
        plt.savefig(fig_path)
    plt.show()
    return

# plot_fig(log_temp = False, save_fig=True)
plot_vel(log_temp = False, save_fig=True)
print('Done')