"""
Code for running an opendrift simulation

The 'opendrift' virtual environment must be activated

There must be an accompanying config.py file in the config_dir where the model is run

Anything in here which we want to make configurable should be moved to the config.py file

"""

import sys, importlib
import numpy as np
from datetime import datetime, timedelta
import opendrift_tools.preprocess as pre_od

# -----------
# OCEANDRIFT
# -----------
def oceandrift(config_dir):
    
    # --------
    # imports
    # --------
    #
    from opendrift.models.oceandrift import OceanDrift
    print('config_dir is '+config_dir)
    sys.path.append(config_dir)
    if 'config' in sys.modules: # this is needed in case we are running this in a loop
        config = importlib.reload(sys.modules['config'])
    else:
        import config
    
    # -------------------------------------
    # initialise openoil and set up readers
    # -------------------------------------
    #
    o = OceanDrift(loglevel=config.loglevel)
    #
    o = pre_od.add_readers(o,config)
    
    # -------------------
    # physical processes
    # -------------------
    #
    # land interaction
    o.set_config('general:use_auto_landmask', True) 
    o.set_config('general:coastline_action', config.coastline_action) 
    o.set_config('general:seafloor_action', config.seafloor_action)
    #
    # prefering to use the exact wind and current...for now anyway
    o.set_config('drift:wind_uncertainty',0) 
    o.set_config('drift:current_uncertainty',0) 
    # ...but add hz diffusivity
    o.set_config('drift:horizontal_diffusivity', config.hz_diff)
    o.set_config('drift:vertical_mixing', config.vert_mix)
    o.set_config('vertical_mixing:timestep', config.vert_mix_tstep)
    o.set_config('drift:vertical_advection', config.vert_adv)
    o.set_config('vertical_mixing:diffusivitymodel','environment') # read from forcing file if it's there
    o.set_config('environment:fallback:ocean_vertical_diffusivity', config.vert_mix_fallback) 
    
    # ------------------
    # seed the elements
    # ------------------
    #
    time_start = datetime.strptime(config.release_start_time, '%Y%m%d_%H')
    time_end = time_start+timedelta(hours=config.release_dur)
    # 
    o.seed_elements(lon=config.lon_release, lat=config.lat_release, z=config.z,
                    radius=config.radius, 
                    time=[time_start,time_end], 
                    number=config.num_part,
                    terminal_velocity = np.random.uniform(config.terminal_min, config.terminal_max, config.num_part),
                    wind_drift_factor=config.wind_drift_factor)
    
    # --------------
    # run the model
    # --------------
    #
    fname = config_dir+'/trajectories.nc' # keeping the filename generic
    o.run(duration=timedelta(days=config.run_dur), time_step=timedelta(minutes=config.time_step), time_step_output=timedelta(minutes=config.time_step_output), outfile=fname) 
    
    # --------
    # cleanup
    # --------
    #
    del(o)
    del(config)
    sys.path.remove(config_dir) # needed if we are running this in a loop



# ----
# OIL
# ----
def oil(config_dir):
    
    # --------
    # imports
    # --------
    #
    from opendrift.models.openoil import OpenOil
    print('config_dir is '+config_dir)
    sys.path.append(config_dir)
    if 'config' in sys.modules: # this is needed in case we are running this in a loop
        config = importlib.reload(sys.modules['config'])
    else:
        import config
    
    # -------------------------------------
    # initialise openoil and set up readers
    # -------------------------------------
    #
    o = OpenOil(loglevel=config.loglevel,  weathering_model='noaa')
    #
    o = pre_od.add_readers(o,config)
    
    # -------------------
    # physical processes
    # -------------------
    #
    # land interaction
    o.set_config('general:use_auto_landmask', True) 
    o.set_config('general:coastline_action', config.coastline_action) 
    o.set_config('general:seafloor_action', config.seafloor_action)
    #
    # prefering to use the exact wind and current...for now anyway
    o.set_config('drift:wind_uncertainty',0) 
    o.set_config('drift:current_uncertainty',0) 
    # ...but add hz diffusivity
    o.set_config('drift:horizontal_diffusivity', config.hz_diff)
    o.set_config('drift:vertical_mixing', config.vert_mix)
    o.set_config('vertical_mixing:timestep', config.vert_mix_tstep)
    o.set_config('drift:vertical_advection', config.vert_adv)
    o.set_config('vertical_mixing:diffusivitymodel','environment') # read from forcing file if it's there
    o.set_config('environment:fallback:ocean_vertical_diffusivity', config.vert_mix_fallback) 
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
    if config.dsd=='uniform':
        o.set_config('seed:droplet_size_distribution','uniform')
        o.set_config('seed:droplet_diameter_min_subsea', config.dsd_param1)
        o.set_config('seed:droplet_diameter_max_subsea', config.dsd_param2)
    elif config.dsd=='lognormal':
        o.set_config('seed:droplet_size_distribution','lognormal')
        o.set_config('seed:droplet_diameter_mu',config.dsd_param1)
        o.set_config('seed:droplet_diameter_sigma',config.dsd_param2)
    else:
        print(config.dsd + 'is not a supported initial droplet size distribution')
    
    time_start = datetime.strptime(config.release_start_time, '%Y%m%d_%H')
    time_end = time_start+timedelta(hours=config.release_dur)
    #
    # flow rate
    o.set_config('seed:m3_per_hour', config.oil_flow_rate); 
    
    # ------------------
    # seed the elements
    # ------------------
    #
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
    
    # --------
    # cleanup
    # --------
    #
    del(o)
    del(config)
    sys.path.remove(config_dir) # needed if we are running this in a loop



# -------
# LEEWAY
# -------
def leeway(config_dir):
    
    # --------
    # imports
    # --------
    #
    from opendrift.models.leeway import Leeway
    print('config_dir is '+config_dir)
    sys.path.append(config_dir)
    if 'config' in sys.modules: # this is needed in case we are running this in a loop
        config = importlib.reload(sys.modules['config'])
    else:
        import config
    
    # -------------------------------------
    # initialise Leeway and set up readers
    # -------------------------------------
    #
    lw = Leeway(loglevel=config.loglevel)
    #
    lw = pre_od.add_readers(lw,config)
    
    # -------------------
    # physical processes
    # -------------------
    #
    # land interaction
    lw.set_config('general:use_auto_landmask', True) 
    lw.set_config('general:coastline_action', config.coastline_action) 
    #
    lw.set_config('drift:horizontal_diffusivity', config.hz_diff)
    
    # ------------------
    # seed the elements
    # -------------------
    #
    time_start = datetime.strptime(config.release_start_time, '%Y%m%d_%H')
    lw.seed_elements(lon=config.lon_release, lat=config.lat_release,
                    radius=config.radius, 
                    time=time_start, 
                    number=config.num_part,
                    object_type=config.object_type)
    
    # --------------
    # run the model
    # --------------
    #
    fname = config_dir+'/trajectories.nc' # keeping the filename generic
    lw.run(duration=timedelta(days=config.run_dur), time_step=timedelta(minutes=config.time_step), time_step_output=timedelta(minutes=config.time_step_output), outfile=fname) 
    
    # --------
    # cleanup
    # --------
    #
    del(lw)
    del(config)
    sys.path.remove(config_dir) # needed if we are running this in a loop
    
