# Configuration file for running an OpenOil simulation 
#
from datetime import datetime, timedelta
#
# ---------------------------------------
# configuration name and run date to use
# ---------------------------------------
#
# give a name for your configuration
# output will be written to a directory with this name
config_name = 'Test_Run'
#
# croco run date
# this is over-written as part of the operational workflow
croco_run_date = 'YYYYMMDD_HH'

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
# start time of spill - use local time (UTC+2)
#spill_start_time=datetime.now() # change to whenever the spill should be 
spill_start_time=datetime(2023,7,25,5,30,0)
#
# duration of the release of oil in hours
release_dur=3
#
# volume of oil spilled in m3
# This is not used directly in the model - it's only used here to get the oil flow rate below
# so you can also specify the 'oil_flow_rate' directly and comment 'oil_volume' if that is convenient 
oil_volume=41.71633
# (35 tonnes, using density of 839 kg/m3 for generic diesel)
#
# oil flow rate in m3/hr
oil_flow_rate=oil_volume/release_dur

# -------------
# forcing files
# -------------
#
# you can just comment files which you don't want to include in the forcing
#
# the reference datetime used in setting up the croco simulations (you shouldn't have to ever change this, unless we reconfigure it in the croco preprocessing)
croco_ref_time = datetime(2000,1,1)
#
# this is an array of file names to allow for the inclusion of multiple croco runs
# The order is important - preference will be given to those which appear first in the array
# The default locations are those insude the docker image used to run operationally 
croco_files = ['/tmp/algoa_01/croco_v1.3.1/C01_I99_MERCATOR_GFS/output/croco_avg.nc',
        '/tmp/swcape_02/croco_v1.3.1/C01_I99_MERCATOR_GFS/output/croco_avg.nc'
        ]
#
# mercator file, as downloaded using the somisana pre-processing tools
mercator_file = '/tmp/downloaded_data/MERCATOR/mercator_YYYYMMDD_HH.nc'

# gfs file, as produced by the croco pre-processing tools after processing the downloaded grb files into a single nc file
gfs_file = '/tmp/downloaded_data/GFS/for_croco/GFS_YYYYMMDD_HH.nc'

# ------------------
# numerical settings
# ------------------
#
# run duration in days
# limit this automatically to avoid run crash at the end of the available model data
run_dur = 10
run_time_max = datetime.strptime(croco_run_date, '%Y%m%d_%H')+timedelta(days=5)
run_dur_max = (run_time_max - spill_start_time).total_seconds()/86400
run_dur=min(run_dur,run_dur_max)
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

# -------------
# grdding info
# -------------
#
# placeholder for now - we can have options for gridding the output here 
# e.g. what grid size to use, what spatial extent

# --------------
# plotting info
# --------------
plot_extents=[25.5,26.5,-34.2,-33.6] # [lon1,lon2,lat1,lat2]
figsize=(8,4) # resize as needed to match the shape of extents below
time_x=0.1 # placement of time label, in axes coordinates
time_y=0.9
vmin=-50   # the z variable is animated so this is your max depth
cmap='Spectral_r' #'Spectral_r' 
plot_cbar=True #True
cbar_loc=(0.9, 0.15, 0.01, 0.7)
croco_dirs_plot=None # croco_dirs
