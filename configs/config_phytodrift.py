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
config_name='phytodrift'
#
# define the date when the croco runs were initialised, in format YYYYMMDD_HH 
run_date='20260507_00'
#
# --------------------------------
# Model Domain
# --------------------------------
#
domain=[14,20,-36,-29]
#
# -------------
# release info
# -------------
# path to phytoplankton flags.
flag_file='/mnt/tmp/downloaded_data/PHYTOPLANKTON_FLAGS/PHYTO_'+run_date+'.nc'           
#
# start time of release, in format YYYYMMDD_HH, in UTC
release_start_time='20260507_00'
#
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
#
# switch to turn on/off use of CROCO as input
use_croco=True 
#
croco_files = ['/mnt/tmp2/sa_southeast_01/croco_v1.3.1/C06_I99_OGCM_WIND_TPXO10/output/croco_avg.nc',
               '/mnt/tmp2/sa_west_02/croco_v1.3.1/C06_I99_OGCM_WIND_TPXO10/output/croco_avg.nc'
               ]
# switch to turn on/off use of OGCM as input
use_ogcm=True
#
# ogcm file, as downloaded using the somisana pre-processing tools
ogcm_files = ['/mnt/tmp2/downloaded_data/OGCM/OGCM_'+run_date+'.nc']
#
# switch to turn on/off use of wind as input
use_wind=False
#
wind_files = ['/mnt/tmp/downloaded_data/WIND/for_croco/U-component_of_wind_Y9999M1.nc',
              '/mnt/tmp/downloaded_data/WIND/for_croco/V-component_of_wind_Y9999M1.nc'
              ]
#
use_waves=True
#
wave_files = ['/mnt/tmp/downloaded_data/CMEMS/CMEMS_'+run_date+'.nc'
              ]
#
# -------------------
# physical processes
# -------------------
#
# coastline interaction - 'stranding' or 'previous'
coastline_action='previous'
#
# seafloor_action - ‘previous’: particles are moved back to previous location 
#                 - ‘deactivate’: particles are deactivated 
#                 - ‘lift_to_seafloor’: particles are lifted vertically to seafloor level
seafloor_action='previous'
#
# include vertical advection
vert_adv=True
#
# include vertical mixing? (boolean)
vert_mix=True
#
# what to use in case the forcing files don't have vertical diffusivity (m2/s)
vert_mix_fallback=0.001
#
# constant horizontal diffusivity (m2/s)
hz_diff = 1
#
# ------------------
# numerical settings
# ------------------
#
# what level of logging to write to the console
# loglevel=50 turns off logs, loglevel=0 gives max information for debugging
loglevel = 0
#
# run duration in days
# default is None, in which case it will be dynamically defined based on the input forcing
# if explicitly defined, make sure the run duration doesn't exceed the temporal range of your inputs!
run_dur = None 
#
# opendrift timestep for particle integration in minutes
# negative value should make the model run in reverse
time_step=15
#
# vertical mixing tstep in seconds
vert_mix_tstep=60
#
# output timestep in minutes
time_step_output=60
