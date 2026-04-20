#!/usr/bin/env python
"""
some generic functions for running opendrift on the operational somisana ocean model output
specifically around setting up the environmental data readers
These are common to all OpenDrift applications

"""

from datetime import datetime
from netCDF4 import num2date
import xarray as xr
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask

def get_seed_points(file_in,domain=None):
    #
    # Function to get the longoitude and latitude positions at points that are flagged 
    # as a red tide in the phytoplankton flag file. 
    #
    if domain is not None:
        ds_flags = ds_flags = xr.open_dataset(file_in).sel(x=slice(domain[0],domain[1]),
                                                           y=slice(domain[3],domain[2])
                                                           )
    else:
        ds_flags = xr.open_dataset(file_in)
            
    phytoplankton_flags = ds_flags.phytoplankton.values.squeeze()
    lons = ds_flags.x.values
    lats = ds_flags.y.values
    LONS,LATS = np.meshgrid(lons,lats)
    ds_flags.close()
    
    idx = np.argwhere((phytoplankton_flags==6))
    j,i = idx[:,0],idx[:,1]
    
    lonp,latp = LONS[j,i],LATS[j,i]
    
    return lonp,latp

def set_croco_time(reader_filename,date_ref):
    # hacky solution to correct the time, as native croco files do not contain reference time
    # code taken directly from reader_ROMS_native.py, just with updated time_units
    print('Setting CROCO time as seconds relative to '+date_ref.strftime('%Y-%m-%d'))
    ocean_time = reader_filename.Dataset.variables['time']
    time_units = 'seconds since '+str(date_ref) 
    reader_filename.times=num2date(ocean_time[:], time_units)
    reader_filename.start_time = reader_filename.times[0]
    reader_filename.end_time = reader_filename.times[-1]
    return reader_filename

def add_readers(o,config):
    # add the somisana readers
    # o - OpenDrift object to add the readers to
    # config - configuration file which contains info on what readers to add
    
    # start by increasing o.max_speed from the default of 1.3 m/s due to warning messages - the Agulhas is fast!
    # we need to make sure we grab enough data from the readers at each time-step (read the docs for what this does)
    # (tests indicate this only works if done before adding the readers)
    o.max_speed = 5 
    
    # -----
    # Land
    # -----
    reader_landmask = reader_global_landmask.Reader()
    o.add_reader(reader_landmask)
    
    # Now start adding the ocean, wind and/or waves input. The order that you do this matters!
    
    # -----------------------------
    # CROCO files covering the run
    # -----------------------------
    #
    # use the reader_ROMS_native reader
    if getattr(config, "use_croco", False):
        print('\nAdding CROCO\n')
        croco_files = config.croco_files
        for croco_file in croco_files:
            reader_croco = reader_ROMS_native.Reader(croco_file)
            croco_ref_time = datetime(config.croco_Yorig,1,1)
            reader_croco = set_croco_time(reader_croco,croco_ref_time)
            o.add_reader(reader_croco)
            print('\n We added CROCO successfully. \n')
    else:
        print('CROCO file(s) not defined - running without CROCO input')
    
    # -------------------------------
    # currents from global OGCM model
    # -------------------------------
    #
    if getattr(config, "use_ogcm", False):
        print('\nAdding OGCM\n')
        ogcm_files = config.ogcm_files
        for ogcm_file in ogcm_files:
            reader_ogcm = reader_netCDF_CF_generic.Reader(ogcm_file)
            o.add_reader(reader_ogcm)
            print('\n We added OGCM successfully. \n')
    else:
        print('OGCM file(s) not defined - running without OGCM input')
        # if you want to exclude for debugging:
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
    
    # -------------
    # Wind forcing
    # -------------
    #
    if getattr(config, "use_wind", False):
        print('\nAdding winds\n')
        wind_files = config.wind_files
        for wind_file in wind_files:
            # Assume we're using the netcdf file on the native grid created during croco preprocessing
            # (but I'm pretty sure this would work on any other regular grid file with standard names)
            # I'm opening the file as an xarray dataset before passing to reader_netCDF_CF_generic
            # because it had trouble with the time conversion inside the reader
            # Doing it like this is a hack to get around this issue, as the time gets handled by xarray
            # rather than by netCDF4's num2date function as done in the reader 
            Dataset = xr.open_mfdataset(wind_file, decode_times=True) # decode_times=True is the default 
            reader_wind = reader_netCDF_CF_generic.Reader(Dataset)    
            o.add_reader(reader_wind)
            print('\n We added winds successfully. \n')
    else:
        print('Wind forcing not defined - running without wind input')
        # if you want to exclude wind for debugging:
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)

    # -------------
    # Wave forcing
    # -------------
    #    
    if getattr(config, "use_waves", False):
        print('\nAdding waves\n')
        waves_files = config.wave_files
        for wave_file in waves_files:
            reader_waves = reader_netCDF_CF_generic.Reader(wave_file)
            o.add_reader(reader_waves)
            print('\n We added waves successfully. \n')
    else:
        print('Wave file(s) not defined - running without wave input')

    return o
