# Configuration file for running an OceanDrift simulation 
#
# We are intentially excluding any python package imports in this configuration file 
#
# The options here are far from exhaustive, but are a few which we've decided to make configurable
# feel free to add more options here, and adjust the run.py file accordingly to read these options
#
# --------------------------------
# configuration name and run date
# --------------------------------
# (this section is only applicable for the operational opendrift runs)
# (these are not used in local runs)
#
# give a name for your configuration
config_name='Test_Run'
#
# define the date when the croco runs were initialised, in format YYYYMMDD_HH 
run_date='20240701_06'

# -------------
# release info
# -------------
#
# coordinates of the release (in geographical degrees)
lon_release=25.74
lat_release=-33.855
#
# depth of the release
# for a subsurface release you can also specify a distance off the seabed like z='seafloor+100' for 100m off the bottom
z=0
#
# radius to be used in initialising the particles
# particles will be initialised around 'lon_release,lat_release' using a standard deviation of 'radius'
# this allows for some initial spreading at location of the initialised particles 
radius=5
#
# start time of spill, in format YYYYMMDD_HH, in UTC
release_start_time='20240701_06'
#
# duration of the release of particles in hours
release_dur=3

# -------------
# forcing files
# -------------
#
# you can just comment files which you don't want to include in the forcing
# the operational workflow uses a sed replcements for OGCM and WIND below to change 
# them to the strings specified in the operational run e.g. MERCATOR and GFS
# For a local run, you of course need to edit the strings to point to the files
# you are forcing with
#
# the Yorig variable used in setting up the croco simulations (used for getting croco file time into real datetimes)
croco_Yorig=2000
#
# this is an array of file names to allow for the inclusion of multiple croco runs
# The order is important - preference will be given to those which appear first in the array
# The default locations are those inside the docker image used to run operationally 
croco_files = ['/mnt/tmp/algoa_01/croco_v1.3.1/C01_I99_OGCM_WIND/output/croco_avg.nc',
        '/mnt/tmp/swcape_02/croco_v1.3.1/C01_I99_OGCM_WIND/output/croco_avg.nc'
        ]

# ogcm file, as downloaded using the somisana pre-processing tools
ogcm_file = '/mnt/tmp/downloaded_data/OGCM/OGCM_'+run_date+'.nc'

# atmospheric forcing file, as produced by the croco pre-processing tools prior to interpolating onto the croco model grid
wind_file = '/mnt/tmp/downloaded_data/WIND/for_croco/WIND_'+run_date+'.nc'

# -------------------
# physical processes
# -------------------
#
# coastline interaction - 'stranding' or 'previous'
coastline_action='stranding'
#
# seafloor_action - ‘previous’: particles are moved back to previous location - ‘deactivate’: particles are deactivated - ‘lift_to_seafloor’: particles are lifted vertically to seafloor level
seafloor_action='lift_to_seafloor'
#
# include vertical mixing? (boolean)
vert_mix=False
#
# what to use in case the forcing files don't have vertical diffusivity (m2/s)
vert_mix_fallback=0.001
#
# inclide vertical advection? (boolean)
vert_adv=False
#
# terminal velocity
# this is an option to specify your own terminal velocity e.g. due to particle negative or positive buoyancy 
# by default we specify a uniform distribution between two values
# negative (positive) values imply sinking (floating) velocities in m/s
terminal_min=0.
terminal_max=0.
#
# constant horizontal diffusivity (m2/s)
hz_diff = 1
#
# wind drift factor
# fraction of the 10 m wind speed used to advect surface particles
wind_drift_factor=0.03

# ------------------
# numerical settings
# ------------------
#
# what level of logging to write to the console
# loglevel=50 turns off logs, loglevel=0 gives max information for debugging
loglevel = 50
#
# run duration in days
# make sure the run duration doesn't exceed the temporal range of your inputs!
run_dur = 0.25
#
# number of particles to release
# generally the more the better, but there are computational limits
# the more particles, the smoother the result will be
# you can calculate the volume of oil per particle upon release as oil_volume/num_part
num_part=5000
#
# opendrift timestep for particle integration in minutes
time_step=15
#
# vertical mixing tstep in seconds
vert_mix_tstep=60
#
# output timestep in minutes
time_step_output=60

# -----------------
# gridding options
# -----------------
#
# options for converting partcle locations into a eulerian grid of particle concentrations
# (used in cli.py for calling the grid_particles function as part of the operational workflow
# for local work, it'll be easier to call the grid_particles function directly in your own python script)
#
grid_type='density' # what kind of gridding to do. Options are 'density', 'surface_oil', 'stranded_oil'
fname_gridded='gridded_density.nc'
grid_extents=None # the spatial extent of the grid [lon0,lon1,lat0,lat1]. If None, then this is automatically determined from the geographic extent of the particles
dx_m=None # grid size in meters, if None, then a 100 x 100 regular grid is generated
max_only=False # option to only write the maximum over the entire file to save disk space (boolean)

# -----------------
# plotting options
# -----------------
#
# options for doing some standardised animations as part of the operational work flow
# the options here are non-exhaustive and are mostly related to the sizing of the plot
# (used in cli.py for calling the plot_particles and plot_gridded functions as part of the operational workflow
# for local work, it'll be easier to call the functions directly in your own python script)
#
fname='trajectories.nc'
# options related to the figure layout
figsize=(8,4) # (hz,vt)
plot_extents = [25.5,26.5,-34.1,-33.6] # spatial extent to plot [lon0,lon1,lat0,lat1]
lscale = 'h' # resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
size_release = 50
# options relating to the dispaly of data
size_scat = 20 # size of the scatter data to be plotted
# by default the particles will just be plotted as black
# options related to the animation
write_gif=True
skip_time = 1 # every nth time-step will be animated (if provided)
tstep_end=None # The last timestep to animate. Only used if write_gif = True. If None, then it'll animate to the end of the file
gif_out_particles='trajectories.gif' # filename of the gif file
gif_out_gridded='gridded_density.gif' # filename of the gif file

