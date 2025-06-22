'''
AUTHOR: Adam Werkema 
DATE CREATED: 5/20/22
DATE LAST EDITED: 6/22/25

DESCRIPTION: This script tracks supercell mesocyclones from CM1 output.
The output from this script provides the location of the mesocyclone at
every output time step, allowing mesocyclone-relative statistics to be
computed. Mesocyclones are tracked by identifying the local maximum of
a Guassian-smoothed w*zeta field at each output time.

INPUT: File must be in netCDF format. This algorithm has not been 
tested on horizontally stretched grids (vertically stretched grids are
fine). This script has been created for simulations with 1-minute 
output from CM1. It has not been tested with output of other frequency.
This script assumes that the domain is in kilometers with (0,0) being
the center.

OUTPUT: Three .npy files that are described below.

meso_key.npy:
Size is [number of output times, number of mesocyclones tracked]. This 
variable is boolean. A '1' means that a mesocyclone is present for the 
given supercell and time. A '0' means that a mesocyclone in not present.

meso_map.npy:
Size is [number of output times, y indices, x indices].
This variable is a 2D map of the mesocyclone(s) at each output time. A 
value of '0' means that the grid cell is not within a mesocyclone. A 
value of '1' means that the grid cell is within mesocyclone #1. A value 
of '2' means that the grid cell is within mesocyclone #2, etc.

meso_centers.npy:
Size is [number of output times, number of mesocyclones tracked, the y 
and x coordinates (in this order)]. This variable tells the estimated 
center of the mesocyclone at each time. Coordinates are given in km
relative to the center of the domain with (0,0) being the center.
'''
###########################################################################
# %% USER SETTINGS

# Select netCDF file
infile = "/home/adamw/data2/multistorm/40_percent_production_runs/30km0deg.nc"
"directory/to/cm1/netCDF/file.nc"

# TIMING VARIABLES

# Time to start tracking mesocyclone(s) (in minutes)
start_time = 40 

# Time to stop tracking mesocyclone (in minutes) (typically the end of the
# simulation)
end_time = 180 

# TUNABLE ALGORITHM VARIABLES 

# Select number of mesocyclones to be tracked
num_storms = 2

# Height to track mesocyclone (in km) (default is 3 km)
algorithm_height = 3

# Minimum w at "algorithm_height" to be considered a mesocyclone (in m/s) 
# (default is 15 m/s)
w_threshold = 15 

# This is the search radius the algorithm uses for finding a mesocyclone at the 
# next output time. This value will be the maximum distance jump the algorithm 
# will make between output times. Make sure this value is large enough so that 
# the algorithm can follow the supercell if the storm is not quasi-stationary 
# in the domain. However, make sure the value isn't too high so that the 
# algorithm will jump to another supercell. (in km) (defaulty is 2.5 km)
algorithm_search_radius = 2.5 

# Fixed radius of the algorithm-identified mesocyclone. This value does not
# affect how the algorithm works, but rather changes the size of the circle
# plotted on the "meso_map". (in km) (default is 5 km)
plotted_meso_radius = 5 

###########################################################################
# %% IMPORT PACKAGES

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib.colors as colors
import matplotlib.animation as animation
import matplotlib.cm as cm
import scipy.ndimage as ndimage

###########################################################################
# %% DEFINE FUNCTIONS 

# Converts from km to an index in python. "ref_array" is "xh", "yh", or "zh".
def km_to_index(km,ref_array):
    index = np.abs(ref_array-km)
    index = np.where(index == np.amin(index))[0][0]
    return index

# Using a radius and grid point, this function plots a circle on a 
# 2D map and returns the map. 'lat_coord', 'lon_coord', and 'radius' are in km. 
# map_2D is 2D array of zeros that has the dimensions of nx and ny.
def make_circle(lat_coord,lon_coord,radius,map_2D):
    
    # Calculate the spatial limits of the circle (in km)
    northernmost = lat_coord + radius
    southernmost = lat_coord - radius
    easternmost = lon_coord + radius
    westernmost = lon_coord - radius
    
    # Convert from metric coordinates (km) to python grid coordinates
    northernmost_index = km_to_index(northernmost,yh) + 1
    southernmost_index = km_to_index(southernmost,yh) - 1
    easternmost_index = km_to_index(easternmost,xh) + 1
    westernmost_index = km_to_index(westernmost,xh) - 1
    
    # Now plot circle to map. Loop through the northernmost/southernmost/etc 
    # indices and select which grid cells will be included in the circle.
    for lat_counter in range(southernmost_index,northernmost_index+1):
        for lon_counter in range(westernmost_index,easternmost_index+1):
            
            # Use distance formula
            distance = ((yh[lat_counter]-lat_coord)**2 + \
                        (xh[lon_counter]-lon_coord)**2)**0.5
            
            # If the grid cell is within the radius of circle, plot it on the 
            # map.
            if distance <= radius:
                map_2D[lat_counter,lon_counter] = True

    # Return the map
    return map_2D.astype(bool)

# Label mesocyclone the 2D "meso_map" if a mesocyclone exists 
# (given the w_threshold criteria specified above). Also, update the meso_key 
# if a mesocyclone exists.
def make_meso_official(uv_max_loc_lat, uv_max_loc_lon, storm_counter):
    
    # Make circle to highlight the area that will be searched for a mesocyclone
    search_circle = make_circle(yh[uv_max_loc_lat], xh[uv_max_loc_lon], 
                                algorithm_search_radius,
                                np.zeros((len(yh),len(xh))))  
    
    # If the maximum value of w within the "plotted_meso_radius" value is 
    # greater than the threshold, document mesocyclone. If not, stop tracking
    # the mesocyclone.
    if np.max(w[time_counter,:,:][search_circle]) >= w_threshold:
        
        # Update the meso_key
        meso_key[time_counter,storm_counter] = 1
        
        # Record the lat/lon coordinates for meso_centers
        meso_centers[time_counter,storm_counter,0] = yh[uv_max_loc_lat]
        meso_centers[time_counter,storm_counter,1] = yh[uv_max_loc_lon] 
        
        # Label on the "meso_map"
        plot_circle = make_circle(yh[uv_max_loc_lat], xh[uv_max_loc_lon], 
                                plotted_meso_radius,
                                np.zeros((len(yh),len(xh)))) 

        meso_map[time_counter,:,:][plot_circle] = storm_counter + 1
        
        # Eliminate values from uv so that the same mesocyclone isn't tracked
        # a second time
        uv[time_counter,:,:][plot_circle] = 0

###########################################################################
# %% LOAD NETCDF DATA

# Point to netCDF file
ds = nc.Dataset(infile)

# Read in coordinate data
xh = ds['xh'][:]
yh = ds['yh'][:]
zh = ds['zh'][:]

# Read in needed fields using height indices
w = ds['winterp'][start_time:end_time+1,km_to_index(algorithm_height, zh),:,:]
zeta = ds['zvort'][start_time:end_time+1,km_to_index(algorithm_height, zh),:,:]

###########################################################################
# %% CALCULATE SMOOTHED UPDRAFT VORTICITY (UV) FOR LAYER

# Allocate 3D array for uv
uv = np.zeros(np.shape(w))

for time_counter in range(end_time-start_time+1):
    
    # Make uv 2D array
    w_slice = ndimage.gaussian_filter(w[time_counter,:,:],10)
    w_slice[w_slice < 0] = 0
    z_slice = ndimage.gaussian_filter(zeta[time_counter,:,:],10)
    z_slice[z_slice < 0] = 0
    uv[time_counter,:,:] =  w_slice * z_slice
    
# Save another copy of uv for plotting
uv_for_plotting = copy.deepcopy(uv)

###########################################################################
# %% MAKE VARIABLES NEEDED FOR TRACKING ALGORITHM
    
# Tells us if a given storm and time has a mesocyclone (1) or not (0)
meso_key = np.zeros((end_time-start_time+1, num_storms)) 

# Map that labels each mesocyclone as 1, 2, 3, etc.
meso_map = np.zeros(np.shape(zeta)) 

# Stores the centroid of the mesocyclone in metric coordinates (for third 
# coordinate, latitude/row comes first)
meso_centers = np.zeros((end_time-start_time+1, num_storms,2)) 

###########################################################################
# %% TRACK MESOCYCLONES

# Loop through each time to determine where mesocyclone(s) are
for time_counter in range(end_time-start_time+1):
    
    # If this is the first time in the simulation, there is no previous 
    # mesocyclone to align with. So, we find the mesocyclone(s).
    if time_counter == 0:
    
        # Initialize the same number of mesocyclones as in "num_storms" 
        for storm_counter in range(num_storms):
            
            # Determine if a uv max in the remaining domain is a mesocyclone. 
            # If so, update the mesocyclone map and key.
            uv_max_loc = np.where(uv[0,:,:] == np.max(uv[0,:,:]))
            make_meso_official(uv_max_loc[0][0], uv_max_loc[1][0], 
                               storm_counter)                                      
        
    # If this isn't the first time in the simulation:    
    else:
        
        # Track the same number of mesocyclones as in "num_storms" 
        for storm_counter in range(num_storms):
            
            # If the mesocyclone has gone extinct at the previous time, skip 
            # search
            if meso_key[time_counter-1,storm_counter] == 0:
                continue
            
            # Print progress to screen
            print("storm number ", storm_counter+1, " at ", 
                  (time_counter+start_time)," minutes")
            
            # Make circle of mesocyclone location indices that is the same 
            # radius as "algorithm_search_radius" around the mesocyclone center
            # for the previous output time
            previous_meso_loc = \
                make_circle(meso_centers[time_counter-1,storm_counter,0],
                            meso_centers[time_counter-1,storm_counter,1],
                            algorithm_search_radius, 
                            np.zeros((len(yh),len(xh))))
    
            # Find the current maximum uv within the previous mesocyclone
            # circle
            uv_max_loc = \
                np.where(uv[time_counter,:,:] == \
                         np.max(uv[time_counter,:,:][previous_meso_loc]))
            
            make_meso_official(uv_max_loc[0][0], uv_max_loc[1][0], 
                               storm_counter)

###########################################################################
# %% SAVE THE OUTPUT ARRAYS

np.save("path/to/output/files/meso_key.npy",meso_key)
np.save("path/to/output/files/meso_map.npy",meso_map)
np.save("path/to/output/files/meso_centers.npy",meso_centers)

###########################################################################
# THE REMAINING CODE IN THIS FILE IS FOR PLOTTING PURPOSES TO CHECK IF
# MESOCYCLONES ARE BEING TRACKED AS INTENDED
###########################################################################       
# %% FIND WESTERNMOST POINT OF UPDRAFTS FOR PLOT LABELING PURPOSES

# Labeler variable for updrafts: 
# [storm number, time, x & y]
updraft_label_location = np.zeros((num_storms, np.shape(w)[0], 2))

# Loop through each storm
for storm_counter in range(num_storms):
    
    # Loop through each time
    for time_counter in range(end_time-start_time+1):
        
        # Make sure mesocyclone exists, otherwise don't assign label
        if meso_key[time_counter,storm_counter] == 1:
            
            # Assign label location
            updraft_label_location[storm_counter, time_counter, 0] = \
                meso_centers[time_counter, storm_counter, 1] - 15
            updraft_label_location[storm_counter, time_counter, 1] = \
                meso_centers[time_counter, storm_counter, 0] - 5
            
        else:
             # Assign x-index
            updraft_label_location[storm_counter,time_counter,0] = 'nan'
            
            # Assign y-index
            updraft_label_location[storm_counter,time_counter,1] = 'nan'

########################################################################### 
# %% PLOT FIGURE FOR SINGLE OUTPUT TIME

# Output time in minutes
time = 120

# Open a figure
fig, axs = plt.subplots(figsize=(10,10))
axs.set_aspect('equal', adjustable='box')
fig.colorbar(cm.ScalarMappable(norm = colors.Normalize(15,75), 
                               cmap = 'gist_ncar'), 
             ax = axs, 
             ticks = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])

# Create plotting coordinates mesh
xx, yy = np.meshgrid(xh,yh)

# Plot reflectivity
reflectivity = ds['dbz'][time, 0,:,:]
axs.contourf(xx, yy, reflectivity, np.arange(20, 76), cmap = 'gist_ncar')

# Label plot
axs.set_xlabel('Zonal Distance (km)')
axs.set_ylabel('Meridional Distance (km)')
axs.set_title("Mesocyclone Tracking Algorithm (blue circle)\n" + \
              "Surface Reflectivity (shaded)\n" + \
              str(algorithm_height) + r" km AGL smoothed w*$\zeta$ " + \
              "(contoured)")

# Make sure mesocyclone exists, otherwise don't plot algorithm
for storm_counter in range(num_storms):

    if meso_key[time-start_time,storm_counter] == 1:
        
        axs.contourf(xx, yy, meso_map[time-start_time,:,:],
                     [storm_counter+0.5,storm_counter+1.5],
                     colors = 'cornflowerblue')
       
        axs.text(updraft_label_location[storm_counter,time-start_time,0],
                 updraft_label_location[storm_counter,time-start_time,1],
                 str(storm_counter+1), size = 20)

# Plot uv contours
if np.max(uv_for_plotting[time-start_time,:,:]) > 0.005:
    axs.contour(xx, yy, uv_for_plotting[time-start_time,:,:],[0.001], 
                colors = 'black')

# Make time label
time = str((time-start_time) + start_time) + " minutes"
axs.text(xh[0]+15, yh[0]+15, time, size = 20)

########################################################################### 
# %% MAKE ANIMATED REFLECTIVITY FIGURE

# Set frames per second
frames_per_sec = 5

# Load in reflectivity data
reflectivity = ds['dbz'][start_time:end_time + 1, 0, :, :]

# Open a figure
fig, axs = plt.subplots(figsize=(10, 10))
axs.set_aspect('equal', adjustable = 'box')
fig.colorbar(cm.ScalarMappable(norm = colors.Normalize(15, 75), 
                               cmap = 'gist_ncar'), 
             ax = axs, 
             ticks = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])

# This is the function that will be called for each frame in the gif
def animate(i):
    
    # Clear plot from previous time
    axs.clear()
    
    # Plot reflectivity
    axs.contourf(xx, yy, reflectivity[i,:,:], np.arange(20, 76), 
                 cmap = 'gist_ncar')
    
    # Label plot
    axs.set_xlabel('Zonal Distance (km)')
    axs.set_ylabel('Meridional Distance (km)')
    axs.set_title("Mesocyclone Tracking Algorithm (blue circle)\n" + \
                  "Surface Reflectivity (shaded)\n" + \
                  str(algorithm_height) + r" km AGL smoothed w*$\zeta$ " + \
                  "(contoured)")
    
    # Make sure mesocyclone exists, otherwise don't plot algorithm
    for storm_counter in range(num_storms):
    
        if meso_key[i, storm_counter] == 1:
            
            axs.contourf(xx, yy, meso_map[i, :, :], 
                         [storm_counter + 0.5, storm_counter + 1.5],
                         colors = 'cornflowerblue')
           
            axs.text(updraft_label_location[storm_counter, i, 0],
                     updraft_label_location[storm_counter, i, 1],
                     str(storm_counter + 1), size = 20)  
            
    # Plot uv contours
    if np.max(uv_for_plotting[i, :, :]) > 0.005:
        axs.contour(xx, yy, uv_for_plotting[i, :, :], [0.001], 
                    colors = 'black')       
            
    # Make time label
    time = str(i + start_time) + " minutes"
    axs.text(xh[0] + 15, yh[0] + 15, time, size = 20)

# This is the driver code for gif creation
ani = animation.FuncAnimation(fig, animate, frames = end_time - start_time + 1,
                              interval = 0, blit = False)

# Save the gif
ani.save("/path/to/output/animation.gif", writer = 'imagemagick', fps = frames_per_sec)