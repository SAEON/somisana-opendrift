# Configuration file for running an OpenOil simulation 
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
config_name='test_oil_blowout'
#
# define the date when the croco runs were initialised, in format YYYYMMDD_HH 
run_date='20241106_12'

# -----------
# spill info
# -----------
#
# coordinates of the release (in geographical degrees)
lon_release=31.779959
lat_release=-30.539622
#
# depth of the release
# for a subsurface release you can also specify a distance off the seabed like z='seafloor+100' for 100m off the bottom
z='seafloor+100'
#
# radius to be used in initialising the particles
# particles will be initialised around 'lon_release,lat_release' using a standard deviation of 'radius'
# this allows for some initial spreading at location of the initialised particles 
# for a subsea blowout this could be hundereds of meters, but a surface spill it will be small, in the order of meters
radius=1000
#
# specify the oil type - important for weathering properties
# Can choose any oil name from https://adios.orr.noaa.gov/oils/
# Or some Norgegain oils from https://opendrift.github.io/oil_types.html
# Or a few other generic oil types added as part of opendrift, such as 'GENERIC INTERMEDIATE FUEL OIL 180'
oil_type='Norman Wells'
#
# start time of spill, in format YYYYMMDD_HH, in UTC
release_start_time='20241102_00'
#
# duration of the release of oil in hours (this can't be zero!)
release_dur=48
#
# volume of oil spilled in m3
# This is not used directly in the model - it's only used here to get the oil flow rate below
# so you can also specify the 'oil_flow_rate' directly and comment 'oil_volume' if that is convenient 
oil_volume=50
#
# oil flow rate in m3/hr
# oil_flow_rate=oil_volume/release_dur
oil_flow_rate=40
#
# initial droplet sizes
# this get's over-ridden by the wave entrainment DSD when oil is submerged into the water column due to wave entrainment
# so this input is not important for a surface release
dsd='lognormal' # 'uniform' or 'lognormal'
# if 'uniform' dsd_param1 is the min diamter (m)
# if 'lognormal' dsd_param1 is the mean diameter (m)
dsd_param1=0.001
# if 'uniform' dsd_param2 is the max diameter (m)
# if 'lognormal' dsd_param2 is the standard deviation (m)
dsd_param2=0.001

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
coastline_action='stranding'
#
# seafloor_action - ‘previous’: particles are moved back to previous location - ‘deactivate’: particles are deactivated - ‘lift_to_seafloor’: particles are lifted vertically to seafloor level
seafloor_action='lift_to_seafloor'
#
# include vertical mixing? (boolean)
vert_mix=True
#
# what to use in case the forcing files don't have vertical diffusivity (m2/s)
vert_mix_fallback=0.001
#
# inclide vertical advection? (boolean)
vert_adv=False
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
# default is None, in which case it will be dynamically defined based on the input forcing
# if explicitly defined, make sure the run duration doesn't exceed the temporal range of your inputs!
run_dur = 7 
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
time_step_output=360

