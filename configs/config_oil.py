# Configuration file for running an OpenOil simulation 
#
# We are intentially excluding any python package imports in this file with the 
# intention of making it easier to transition to a gui interface at some stage
#
# --------------------------------
# configuration name and run date
# --------------------------------
#
# give a name for your configuration
# output will be written to a directory with this name
config_name='Test_Run'
#
# define the date when the croco runs were initialised, in format YYYYMMDD_HH 
# (this is only applicable for the operational opendrift runs)
run_date='20240627_06'

# -----------
# spill info
# -----------
#
# coordinates of the spill (in geographical degrees)
lon_spill=25.74
lat_spill=-33.855
#
# depth of the release
# For a surface release I prefer to put a small negative number like z=-0.001
# this effectively means the weathering is applied after mixing in the first time-step (not a big deal in the end)
# for a subsurface release you can also specify a distance off the seabed like z='seafloor+100' for 100m off the bottom
z=-0.001
#
# radius to be used in initialising the particles
# particles will be initialised around 'lon_spill,lat_spill' using a standard deviation of 'radius'
# this allows for some initial spreading at location of the initialised particles 
# for a subsea blowout this could be hundereds of meters, but a surface spill it will be small, in the order of meters
radius=5
#
# specify the oil type - important for weathering properties
# Can choose any oil name from https://adios.orr.noaa.gov/oils/
# Or some Norgegain oils from https://opendrift.github.io/oil_types.html
# Or a few other generic oil types added as part of opendrift, such as 'GENERIC INTERMEDIATE FUEL OIL 180'
oil_type='GENERIC INTERMEDIATE FUEL OIL 180'
#
# start time of spill, in format YYYYMMDD_HH, in UTC
spill_start_time='20240627_06'
#
# duration of the release of oil in hours
release_dur=3
#
# volume of oil spilled in m3
# This is not used directly in the model - it's only used here to get the oil flow rate below
# so you can also specify the 'oil_flow_rate' directly and comment 'oil_volume' if that is convenient 
oil_volume=50
#
# oil flow rate in m3/hr
oil_flow_rate=oil_volume/release_dur

# -------------
# forcing files
# -------------
#
# you can just comment files which you don't want to include in the forcing
#
# the Yorig variable used in setting up the croco simulations (used for getting croco file time into real datetimes)
croco_Yorig=2000
#
# this is an array of file names to allow for the inclusion of multiple croco runs
# The order is important - preference will be given to those which appear first in the array
# The default locations are those insude the docker image used to run operationally 
croco_files = ['/mnt/tmp/algoa_01/croco_v1.3.1/C01_I99_OGCM_WIND/output/croco_avg.nc',
        '/mnt/tmp/swcape_02/croco_v1.3.1/C01_I99_OGCM_WIND/output/croco_avg.nc'
        ]

# ogcm file, as downloaded using the somisana pre-processing tools
# the operational workflow uses sed to replace OGCM with the actual string e.g. MERCATOR
ogcm_file = '/mnt/tmp/downloaded_data/OGCM/OGCM_'+run_date+'.nc'

# atmospheric forcing file, as produced by the croco pre-processing tools prior to interpolating onto the croco model grid
# the operational workflow uses sed to replace WIND with the actual string e.g. GFS
wind_file = '/mnt/tmp/downloaded_data/WIND/for_croco/WIND_'+run_date+'.nc'

# ------------------
# numerical settings
# ------------------
#
# run duration in days
# make sure the run duration doesn't exceed the temporal range of your inputs!
run_dur = 5
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
# output timestep in minutes
time_step_output=60
#
# constant horizontal diffusivity (m2/s)
# this is applied to the random walk component of the particle motion
hz_diff = 1
#
# wind drift factor
# fraction of the 10 m wind speed used to advect surface particles
wind_drift_factor=0.03

# -----------------
# plotting options
# -----------------
#
# options for doing plots as part of the operational work flow
#
fname='/mnt/tmp/opendrift_oil/'+config_name+'/OGCM_WIND/trajectories.nc'
tstep=0 # the step to plot, or the first step to animate.
# options related to the figure layout
figsize=(8,4) # (hz,vt)
extents = [25.5,26.5,-34.1,-33.6] # spatial extent to plot [lon0,lon1,lat0,lat1]
lscale = 'h' # resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
# options relating to the release location
lon_release=lon_spill
lat_release=lat_spill
size_release = 50
# options relating to the dispaly of data, colormap and colorbar
size_scat = 20 # size of the scatter data to be plotted
var_str='z' # variable to plot
ticks = [-1,0] #[0,5,10,15,20,25,30] # the ticks to plot relating to the colormap (can be irregularly spaced)
cmap = 'Greys' #'Spectral_r' # colormap to use
plot_cbar = False
cbar_loc = [0.85, 0.2, 0.02, 0.6] # where on the plot to put the colorbar [x,y,width,height]
cbar_label = 'depth (m)'
# options related to the plot output file
jpg_out=None # filename of the jpg file
write_jpg=False
# options related to the animation
gif_out='/mnt/tmp/opendrift_oil/'+config_name+'/OGCM_WIND/trajectories_z.gif' # filename of the gif file
write_gif=True
skip_time = 1 # every nth time-step will be animated (if provided)
tstep_end=None # The last timestep to animate. Only used if write_gif = True. If None, then it'll animate to the end of the file


