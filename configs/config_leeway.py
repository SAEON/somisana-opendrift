# Configuration file for running a Leeway simulation 
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
config_name='test_leeway'
#
# define the date when the croco runs were initialised, in format YYYYMMDD_HH 
run_date='20241106_12'

# ------------
# object info
# ------------
#
# coordinates of the release (in geographical degrees)
lon_release=18
lat_release=-34
#
# radius to be used in initialising the particles
# particles will be initialised around 'lon_release,lat_release' using a standard deviation of 'radius'
# this allows for some initial spreading at location of the initialised particles
radius=3000
#
# specify the object type - important for drift properties
# see references here - https://opendrift.github.io/autoapi/opendrift/models/leeway/index.html#module-opendrift.models.leeway
# and object types here - https://github.com/OpenDrift/opendrift/blob/master/opendrift/models/OBJECTPROP.DAT
# user must specify an integer corresponding to the object type number in the table
object_type=26

# start time of spill, in format YYYYMMDD_HH, in UTC
release_start_time='20241102_00'

# -------------
# forcing files
# -------------
#
# you can just comment files which you don't want to include in the forcing
# the operational workflow uses sed replcements for OGCM and WIND below to change 
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
croco_files = ['/mnt/tmp/sa_southeast_01/croco_v1.3.1/C04_I99_OGCM_WIND/output/croco_avg.nc',
        '/mnt/tmp/sa_west_02/croco_v1.3.1/C04_I99_OGCM_WIND/output/croco_avg.nc'
        ]

# ogcm file, as downloaded using the somisana pre-processing tools
ogcm_file = '/mnt/tmp/downloaded_data/OGCM/OGCM_'+run_date+'.nc'

# atmospheric forcing file, as produced by the croco pre-processing tools
wind_files = ['/mnt/tmp/downloaded_data/WIND/for_croco/U-component_of_wind_Y9999M1.nc',
        '/mnt/tmp/downloaded_data/WIND/for_croco/V-component_of_wind_Y9999M1.nc'
        ]

# -------------------
# physical processes
# -------------------
#
# coastline interaction - 'stranding' or 'previous'
coastline_action='previous'
#
# constant horizontal diffusivity (m2/s)
hz_diff = 0

# ------------------
# numerical settings
# ------------------
#
# what level of logging to write to the console
# loglevel=50 turns off logs, loglevel=0 gives max information for debugging
loglevel = 50
#
# run duration in days
# default is None, in which case it will be dynamically defined based on the input forcing
# if explicitly defined, make sure the run duration doesn't exceed the temporal range of your inputs!
run_dur = None 
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

