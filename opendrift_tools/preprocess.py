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
        wind_file=config.wind_file # this will fail if wind_file is not defined
        if 'GFS' in wind_file:
            # Assume we're using the netcdf file on the gfs native grid created during croco preprocessing 
            # I'm opening the file as an xarray dataset before passing to reader_netCDF_CF_generic
            # because it had trouble with the time conversion inside the reader
            # Doing it like this is a hack to get around this issue, as the time gets handled by xarray
            # rather than by netCDF4's num2date function as done in the reader 
            Dataset = xr.open_dataset(wind_file, decode_times=True) # decode_times=True is the default 
            reader_gfs = reader_netCDF_CF_generic.Reader(Dataset,
                                standard_name_mapping={'uwnd': 'x_wind',
                                                       'vwnd': 'y_wind'},)    
            o.add_reader(reader_gfs)
        # add elif statements here for other wind forcings as needed
        # the else statement catches other files which can be read with the generic reader
        else:
            print('no wind forcing defined')
            reader_wind = reader_netCDF_CF_generic.Reader(wind_file)
            o.add_reader(reader_wind)
    except:
        print('no wind forcing defined, running without wind')
        # if you want to exclude wind for debugging:
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)

    return o

def WASA3_2_nc(wasa_grid, 
                 wasa_zarr_dir, 
                 out_dir,
                 start_date,
                 end_date
                 ):
    '''
    Create monthly nc files from WASA3 u10, v10 data
    To be used as OpenDrift forcing
    
    Parameters
    ----------
    wasa_grid     : /your_dir/geo_em.d03.nc
    wasa_zarr_dir : /your_dir/uv_ds_10m
    out_wasa_dir  : the directory to save blk files with the WASA data
    start_date    : start datetime for processing (to the nearest month)
    end_date      : end datetime for processing (to the nearest month)
    
    Returns
    -------
    Writes monthly netcdf files in out_dir

    '''
        
    # get the relevant wasa data
    ds_wasa=xr.open_zarr(wasa_zarr_dir)
    ds_wasa_grd=xr.open_dataset(wasa_grid)
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze()
    
    # extract the wasa grid variables we need from this spatial subset
    lon_wasa = ds_wasa_grd.XLONG_M.squeeze()
    lat_wasa = ds_wasa_grd.XLAT_M.squeeze() 
    cosa = ds_wasa_grd.COSALPHA.squeeze()
    sina = ds_wasa_grd.SINALPHA.squeeze()
    # flatten lon,lat arays for use in interpolation later
    source_points = np.column_stack((np.array(lon_wasa).flatten(), np.array(lat_wasa).flatten()))

    # loop through months
    date_now=start_date
    while date_now <= end_date:
        file_out = 'WASA3'+date_now.strftime('_Y%YM%m')+'.nc'
   
        if os.path.isfile(out_dir+'/'+file_out):
            print(file_out + ' already exists - skipping')
        else:
            print('\nworking on '+file_out)
            
            # start and end days of this month
            start_extract=datetime(date_now.year,date_now.month,1,0,0,0)
            day_end = calendar.monthrange(date_now.year, date_now.month)[1]
            end_extract=datetime(date_now.year,date_now.month,day_end,23,59,59)
            
            # subset the WASA data for this month
            ds_wasa_now = ds_wasa.sel(Time=slice(start_extract,end_extract))
            time_now = ds_wasa_now.Time
            
            # rotate the wasa vectors to be east,north components
            u = ds_wasa_now.U10
            v = ds_wasa_now.V10
            u_wasa = u * cosa - v * sina
            v_wasa = v * cosa + u * sina
            
            for t in range(len(time_now)):
                # Extract wasa values at this specific time-step
                u_wasa_now = u_wasa.isel(Time=t)
                v_wasa_now = v_wasa.isel(Time=t)
                
                # check if there are any missing data and print a warning
                # This is needed due to some random blocks of missing values we've found
                # I'm automatically interpolating over these data so we have workable files
                # THIS NEEDS TO BE SORTED OUT WITH CSAG!
                # Then we can remove this part of the code
                if np.any(np.isnan(u_wasa_now.values.flatten())):
                    print('wasa u10 contains missing data - '+time_now[t])
                    u_wasa_interp = u_wasa_now.copy()
                if np.any(np.isnan(v_wasa_now.values.flatten())):
                    print('wasa v10 contains missing data - '+time_now[t])
            
            # # interpolate the wasa u,v wind vectors onto the croco grid   
            # # I tried to find an in-built xarray interpolation method to do the job
            # # efficiently, but it only works for regular grids
            # # So we've got to extract the data and use scipy's griddata function
            # # Note I'm puposefully interpolating onto the croco rho grid (not the croco u,v grids)
            # # because we have to rotate the u,v vectors (on the rho grid)
            # # before writing the data to file, so would have to use u2rho() and v2rho() 
            # # which would introduce additional unnecessary interpolation
            # u_wasa_interp = np.zeros_like(ds_blk.wspd) # using 'wspd' to get the rho grid
            # v_wasa_interp = np.zeros_like(ds_blk.wspd)
            # for t in range(len(blk_time)):
            #     if t % 50 == 0:
            #         percentage_complete = (t/len(blk_time)) * 100
            #         print(f"{percentage_complete:.0f}% complete")
                
            #     # Extract wasa values at this specific time-step
            #     # flattening for use in griddata
            #     u_wasa_now = u_wasa.isel(Time=t).values.flatten()
            #     v_wasa_now = v_wasa.isel(Time=t).values.flatten()
                
            #     # Perform interpolation 
            #     # 'nearest' method speeds it up a bit
            #     # but I guess 'linear' is better
            #     u_wasa_now_interp = griddata(source_points, u_wasa_now, (lon_rho, lat_rho), method=interp_method) 
            #     v_wasa_now_interp = griddata(source_points, v_wasa_now, (lon_rho, lat_rho), method=interp_method)             
                
            #     # hmmm, with 'linear' interpolation we're left with some nan's in the case where our domain 
            #     # extends a little past the WASA data, so default to nearest wasa data here
            #     # this will slow things down a bit, but not sure what else to do, other than regenerate the croco grid!
            #     if np.any(np.isnan(u_wasa_now_interp)):
            #         u_wasa_now_interp_nearest = griddata(source_points, u_wasa_now, (lon_rho, lat_rho), method='nearest') 
            #         v_wasa_now_interp_nearest = griddata(source_points, v_wasa_now, (lon_rho, lat_rho), method='nearest')
            #         # now replace missing data with nearest available data
            #         u_wasa_now_interp[np.isnan(u_wasa_now_interp)] = u_wasa_now_interp_nearest[np.isnan(u_wasa_now_interp)]
            #         v_wasa_now_interp[np.isnan(v_wasa_now_interp)] = v_wasa_now_interp_nearest[np.isnan(v_wasa_now_interp)]
                
            #     u_wasa_interp[t,:,:] = u_wasa_now_interp
            #     v_wasa_interp[t,:,:] = v_wasa_now_interp
                
            # # calculate the wspd from the wasa components (on the rho grid)
            # spd_wasa_interp = np.hypot(u_wasa_interp, v_wasa_interp)
            
            # # uwnd and vwnd are grid aligned, so we need to rotate vectors
            # # (based on official croco_tools interp_ERA5.m)
            # cos_a = np.cos(angle.values)
            # sin_a = np.sin(angle.values)
            # # have to use the data on the rho grid for the rotation calc 
            # u_out = u_wasa_interp*cos_a + v_wasa_interp*sin_a
            # v_out = v_wasa_interp*cos_a - u_wasa_interp*sin_a
            # # but get onto the u,v grids to write to file
            # u_out = post.rho2u(u_out)
            # v_out = post.rho2v(v_out)
            
            # # Update uwnd in ds_blk
            # ds_blk['uwnd'] = (('bulk_time', 'eta_u', 'xi_u'), u_out)
            # ds_blk['vwnd'] = (('bulk_time', 'eta_v', 'xi_v'), v_out)
            # ds_blk['wspd'] = (('bulk_time', 'eta_rho', 'xi_rho'), spd_wasa_interp)
            
            # # I'm pretty sure we don't need sustr and svstr in blk files as this is computed online
            # # So I'm removing these variables to be sure
            # # leaving them hanging around seems a little messy since we're replacing uwnd and vwnd
            # # while we're at it, radlw is also not used (radlw_in is), so remove that as well
            # ds_blk.drop_vars(["sustr", "svstr", "radlw"])
            
            # # write float32 instead of float64 in an attempt to save space
            # encoding = {
            #     "tair": {"dtype": "float32"},
            #     "rhum": {"dtype": "float32"},
            #     "prate": {"dtype": "float32"},
            #     "wspd": {"dtype": "float32"},
            #     "radlw_in": {"dtype": "float32"},
            #     "radsw": {"dtype": "float32"},
            #     "uwnd": {"dtype": "float32"},
            #     "vwnd": {"dtype": "float32"},
            #     "bulk_time": {"dtype": "float32"},
            # }
            
            # # write the new blk file
            # ds_blk.to_netcdf(out_wasa_dir+'/'+file_out,
            #                  encoding=encoding)
            
            # ds_blk.close()
            
        date_now=date_now+timedelta(days=32) # 32 days ensures we get to the next month
        date_now=datetime(date_now.year, date_now.month, 1) # set the first day of the month 
        
    ds_wasa.close()
    ds_wasa_grd.close()