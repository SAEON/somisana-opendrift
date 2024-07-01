#!/usr/bin/env python
"""
Code for running an opendrift simulation

The 'opendrift' virtual environment must be activated

There must be an accompanying config.py file in the config_dir where the model is run

Anything in here which we want to make configurable should be moved to the config.py file

"""

import sys
from datetime import datetime, timedelta
import opendrift_tools.preprocess as pre_od

def oil(config_dir):
    
    from opendrift.models.openoil import OpenOil
    print('config_dir is '+config_dir)
    sys.path.append(config_dir)
    import config
    # -------------------------------------
    # initialise openoil and set up readers
    # -------------------------------------
    #
    o = OpenOil(loglevel=0,  weathering_model='noaa')  # loglevel=50 turns off logs, loglevel=0 gives max information for debugging
    
    o = pre_od.add_readers(o,config)
    
    # -------------------
    # physical processes
    # -------------------
    #
    # land interaction
    o.set_config('general:use_auto_landmask', True) 
    o.set_config('general:coastline_action', 'stranding') 
    o.set_config('general:seafloor_action', 'lift_to_seafloor')
    #
    # I'd prefer to use the exact wind and current...for now anyway
    o.set_config('drift:wind_uncertainty',0) 
    o.set_config('drift:current_uncertainty',0) 
    # ...but add hz diffusivity
    o.set_config('drift:horizontal_diffusivity', config.hz_diff) # could test sensitivity 
    o.set_config('drift:vertical_mixing', True) 
    o.set_config('vertical_mixing:timestep', 60) # seconds, higher timestep for this process
    o.set_config('drift:vertical_advection', False) # maybe something to test sensitivity to
    o.set_config('vertical_mixing:diffusivitymodel','environment') # default - reading from croco file
    o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.005) # loosely based on near-surface CROCO output. We'll need this for when we combine our croco runs which include diffusivity with mercator data which does not
    o.set_config('wave_entrainment:entrainment_rate', 'Li et al. (2017)') # default
    o.set_config('wave_entrainment:droplet_size_distribution', 'Li et al. (2017)') # as per Rohrs et al. (2018)
    
    # ---------------------
    # oil spill properties
    # ---------------------
    #
    # weathering
    o.set_config('processes:dispersion', False) # not needed as we are modelling wave entrainment explicitly (see Rohrs et al. (2018))
    o.set_config('processes:evaporation', True)
    o.set_config('processes:emulsification', True)
    #
    # droplet sizes
    # this is only important for a subsea release, and I assume we're looking at surface releases for now
    # so I'm just specifying large initial droplets so that all particles can surface after 1 time-step
    # (this get's over-ridden by the wave entrainment DSD when pushed into the water column due to wave entrainment)
    o.set_config('seed:droplet_size_distribution','uniform')
    o.set_config('seed:droplet_diameter_min_subsea', 0.05)
    o.set_config('seed:droplet_diameter_max_subsea', 0.05)
    # if it's a subsea blowout you may want to specify a lognormal distribution
    #o.set_config('seed:droplet_size_distribution','lognormal')
    #o.set_config('seed:droplet_diameter_mu',1.5e-3)
    #o.set_config('seed:droplet_diameter_sigma',0.5e-3)
    #
    time_start = datetime.strptime(config.spill_start_time, '%Y%m%d_%H')
    time_end = time_start+timedelta(hours=config.release_dur)
    #
    # flow rate
    o.set_config('seed:m3_per_hour', config.oil_flow_rate); 
    #
    # seed the elements
    o.seed_elements(lon=config.lon_release, lat=config.lat_release, z=config.z,
                    radius=config.radius, 
                    time=[time_start,time_end], 
                    number=config.num_part,
                    oil_type=config.oil_type,
                    wind_drift_factor=config.wind_drift_factor)
    
    # --------------
    # run the model
    # --------------
    #
    fname = config_dir+'/trajectories.nc' # keeping the filename generic
    o.run(duration=timedelta(days=config.run_dur), time_step=timedelta(minutes=config.time_step), time_step_output=timedelta(minutes=config.time_step_output), outfile=fname) 
     
    del(o)
    
