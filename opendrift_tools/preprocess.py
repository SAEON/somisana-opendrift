#!/usr/bin/env python
"""
some generic functions for running opendrift on the operational somisana ocean model output
specifically around setting up the environmental data readers
These are common to all OpenDrift applications

"""

from netCDF4 import num2date
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask

def set_croco_time(reader_filename,date_ref):
    # hacky solution to correct the time, as native croco files do not contain reference time
    # code taken directly from reader_ROMS_native.py, just with updated time_units
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
    # if you want to exclude currents for debugging:
    #o.set_config('environment:fallback:x_sea_water_velocity', 0)
    #o.set_config('environment:fallback:y_sea_water_velocity', 0)
    #
    # use the reader_ROMS_native reader
    if config.croco_files:
        for croco_file in config.croco_files:
            reader_croco = reader_ROMS_native.Reader(croco_file)
            reader_croco = set_croco_time(reader_croco,config.croco_ref_time)
            o.add_reader(reader_croco)
    
    # ------------------
    # MERCATOR currents
    # ------------------
    #
    if config.mercator_file:
        reader_mercator = reader_netCDF_CF_generic.Reader(config.mercator_file)
        o.add_reader(reader_mercator)
    
    # ---------
    # GFS wind
    # ---------
    #
    # if you want to exclude wind for debugging:
    # o.set_config('environment:fallback:x_wind', 0)
    # o.set_config('environment:fallback:y_wind', 0)
    #
    # I'm opening the file as an xarray dataset before passing to reader_netCDF_CF_generic
    # because it had trouble with the time conversion inside the reader
    # Doing it like this is a hack to get around this issue, as the time gets handled by xarray
    # rather than by netCDF4's num2date function as done in the reader 
    if config.gfs_file:
        Dataset = xr.open_dataset(config.gfs_file, decode_times=True) # decode_times=True is the default 
        reader_gfs = reader_netCDF_CF_generic.Reader(Dataset,
                                standard_name_mapping={'uwnd': 'x_wind',
                                                       'vwnd': 'y_wind'},)    
    o.add_reader(reader_gfs)
    
    return o
