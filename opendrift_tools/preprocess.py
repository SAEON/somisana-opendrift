#!/usr/bin/env python
"""
some generic functions for running opendrift on the operational somisana ocean model output
specifically around setting up the environmental data readers
These are common to all OpenDrift applications

"""

import os
from datetime import datetime, timedelta
import calendar
from netCDF4 import num2date
import xarray as xr
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask

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
    
    # Now start adding the ocean and wind input. The order that you do this matters!
    # Preference is given to data which is added first. We could include boolean checks on 
    # whether certain forcings should be included e.g. SAWS winds, and different ocean forcings
    
    # -----------------------------
    # CROCO files covering the run
    # -----------------------------
    #
    # use the reader_ROMS_native reader
    try:
        croco_files = config.croco_files
        for croco_file in croco_files:
            reader_croco = reader_ROMS_native.Reader(croco_file)
            croco_ref_time = datetime(config.croco_Yorig,1,1)
            reader_croco = set_croco_time(reader_croco,croco_ref_time)
            o.add_reader(reader_croco)
    except:
        print('CROCO currents not defined, running without them')
    
    # -------------------------------
    # currents from global OGCM model
    # -------------------------------
    #
    try:
        ogcm_file=config.ogcm_file
        reader_ogcm = reader_netCDF_CF_generic.Reader(ogcm_file)
        o.add_reader(reader_ogcm)
    except:
        print('ogcm currents not defined, using a fallback of ZERO')
        # if you want to exclude for debugging:
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
    
    # -------------
    # Wind forcing
    # -------------
    #
    try:
        wind_files = config.wind_files
        for wind_file in wind_files:
            reader_wind = reader_netCDF_CF_generic.Reader(wind_file)
            o.add_reader(reader_wind)
    except:
        print('no wind forcing defined, running without wind')
        # if you want to exclude wind for debugging:
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)

    return o
