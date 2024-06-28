"""
some functions for processing raw opendrift output

"""

import os
import numpy as np
import xarray as xr
import opendrift
import xarray as xr
from xhistogram.xarray import histogram

def fill_deactivated(ds):
    '''
    data for deactivated particles are assigned to being missing, so this function
    fills in the missing values at the end of the file with the last valid time-step
    This allows us to more easily plot and analyse the data
    '''
    # use valid longitude values to mask deactivated particles
    ds_masked = ds.where(ds.lon < 361.)
    # use ffill to fill in data using the last valid time-step
    # this won't fill missing data at the start of the file if a particle isn't yet activated
    # but will fill once it has been deactivated
    ds_filled=ds_masked.ffill('time')
    return ds_filled

def get_lonlat_bins(extents, dx_m):
    '''
    create a regular grid from the input extents and a grid size dx_m
    code copied from opendrift/models/basemodel.py 
    '''
    lonmin = extents[0]
    lonmax = extents[1]
    latmin = extents[2]
    latmax = extents[3]
    deltalat = dx_m / 111000.0  # m to degrees
    deltalon = deltalat / np.cos(np.radians((latmin + latmax) / 2))
    latbin = np.arange(latmin - deltalat, latmax + deltalat, deltalat)
    lonbin = np.arange(lonmin - deltalon, lonmax + deltalon, deltalon)
    return lonbin, latbin

def get_stranded_flag(ds):
    '''
    Identify what the flag index is in the file (it is not the same for all files)
    (maybe generalise this function to handle any flag string)
    '''
    # identify stranded particles
    flag_meanings=ds.status.flag_meanings
    if not 'stranded' in flag_meanings:
        # if no stranded particles, add this flag manually - 
        # just so we can use this to do subsetting on stranded particles later 
        flag_meanings=flag_meanings+' stranded'
    stranded_flag=flag_meanings.split().index('stranded')
    return stranded_flag

def subset_tstep(ds,tstep):
    '''
    extract subset of the dataset for a particular timestep
    The only reason for making a dedicated function is so that
    tstep can be either an integer or a datetime object
    '''
    if isinstance(tstep,int):
        time_plot = ds.time.values[tstep]
        ds_time = ds.sel(time=time_plot)
    else:
        ds_time = ds.sel(time=tstep)
    
    return ds_time

def grid_surface_oil(iteration_dir,dx_m=25000,extents=[16,34,-40,-28],time_factor=1):
    
    os.chdir(iteration_dir)
    lonbin, latbin = get_lonlat_bins(extents,dx_m)

    outfile = 'trajectories.nc'
    # use xarray to grid the particles onto a regular grid
    oa = opendrift.open_xarray(outfile)
    ds=oa.ds
    
    # identify stranded particles
    stranded_flag = get_stranded_flag(ds)

    # get data for computing surface thickness
    #
    # maybe use a threshold of say -0.1 m?
    ds_surf=ds.where(ds.z==0.0) # could use drop=True to reduce the size of the data, but should make no difference to end result
    # remove stranded data
    ds_surf=ds_surf.where(~(ds_surf.status==stranded_flag)) 
    #
    # compute the weights input to the histogram function which will tell us 
    # the oil volume (m3) per surface particle
    # it will be multipled by the histogram (which computes the number of particles per grid cell)
    # therefore providing us with oil volume per grid cell
    # note we are using mass of the oil emulsion into account (see openoil.get_oil_budget)
    # mass_emulsion = mass_oil / (1 - water_fraction)
    surf_volume=(ds_surf.mass_oil / (1 - ds_surf.water_fraction))/ds_surf.density  
    #
    # compute surface thickness histogram
    h_surf = histogram(ds_surf.lon,
                  ds_surf.lat,
                  bins=[lonbin, latbin],
                  dim=['trajectory'],
                  weights=surf_volume,
                  density=False)
    # convert oil volume per grid cell into oil thickness in micron
    h_surf = h_surf/(np.power(dx_m,2))*1e6
    #
    # TODO: add minimum time to oiling like we do for the stranding?
    #
    # compute the maximum surface thickness over the run
    # TODO: MAYBE WE WANT TO KEEP THE SURFACE THICKNESS AT EACH TIME-STEP? same applies for stranded
    h_surf_max=h_surf.max(('time'))
    h_surf_max=h_surf_max.rename('maximum')
    h_surf_max.attrs["standard_name"] = 'maximum_surface_oil_thickness'
    h_surf_max.attrs["units"] = 'micron'
    # add the bin edges of the histogram so we can use pcolormesh plot properly
    # just adding it as an attribute
    h_surf_max.attrs["lon_bin_edges"] = lonbin
    h_surf_max.attrs["lat_bin_edges"] = latbin
    h_surf_max.to_netcdf('surface_dx'+str(dx_m)+'m.nc')

def grid_stranded_oil(iteration_dir,dx_m=25000,extents=[16,34,-40,-28],time_factor=1):
    
    os.chdir(iteration_dir)
    lonbin, latbin = get_lonlat_bins(extents,dx_m)

    outfile = 'trajectories.nc'
    # use xarray to grid the particles onto a regular grid
    oa = opendrift.open_xarray(outfile)
    ds=oa.ds
    
    # fill deactivated particles with data from the last active time-step
    # this will allow stranded particles to accumulate instead of being deactivated after stranding
    ds=fill_deactivated(ds)
    
    # identify stranded particles
    stranded_flag = get_stranded_flag(ds)
    
    # compute the stranded concentration
    ds_strand=ds.where(ds.status==stranded_flag)
    # stranded volume in the same way as surface volume
    strand_mass=(ds_strand.mass_oil / (1 - ds_strand.water_fraction)) 
    h_strand = histogram(ds_strand.lon,
                  ds_strand.lat,
                  bins=[lonbin, latbin],
                  dim=['trajectory'],
                  weights=strand_mass,
                  density=False)
    # convert oil volume per grid cell into g/m2
    beach_width_m=30 # m 1.5 m tidal range and 1:20 beach slope
    h_strand = h_strand*1000/dx_m/beach_width_m
    # compute the minimum time to stranding
    days_since_start = (h_strand.time.data-h_strand.time.data[0]).astype('timedelta64[s]').astype(np.int32)/3600/24 # seriously convoluted but it works. Timedelta64 is not ideals
    h_strand_time=xr.full_like(h_strand, 10e6) # using 10e6 as an arbitrary large number
    # loop through time and replace stranded grid cells with time in days since spill start
    for ii, days in enumerate(days_since_start):
        # h_strand_time.data[ii,:,:] = np.where(np.array(h_strand.data[ii,:,:])>0,days,10e6)
        h_strand_time.data[ii,:,:]=xr.where(h_strand.data[ii,:,:]>0,days,10e6)
    # compute the minimum stranding time
    h_strand_time_min=h_strand_time.min(('time'))
    h_strand_time_min=h_strand_time_min.rename('minimum_time')
    h_strand_time_min.attrs["standard_name"] = 'minimum_time_to_stranding'
    h_strand_time_min.attrs["units"] = '_days_' # using underscores to avoid xarray reading this variable as a timedelta64 later
    # compute the maximum stranded concentration over the run
    h_strand_max=h_strand.max(('time'))
    h_strand_max=h_strand_max.rename('maximum')
    h_strand_max.attrs["standard_name"] = 'maximum_stranded_oil_density'
    h_strand_max.attrs["units"] = 'g/m2'
    # add the bin edges of the histogram so we can use pcolormesh plot properly
    # just adding it as an attribute
    h_strand_max.attrs["lon_bin_edges"] = lonbin
    h_strand_max.attrs["lat_bin_edges"] = latbin
    #
    h_strand_merge=xr.merge([h_strand_max, h_strand_time_min])
        
    h_strand_merge.to_netcdf('stranded_dx'+str(dx_m)+'m.nc')

def get_trajectories_oil_budget(iteration_dir):
    os.chdir(iteration_dir)
    outfile = 'trajectories.nc'
    
    # # read the data
    oa = opendrift.open_xarray(outfile)
    ds=oa.ds
    
    # fill deactivated particles with data from the last active time-step (NB for mass budgets)
    ds=fill_deactivated(ds)
    
    # identify stranded particles
    stranded_flag = get_stranded_flag(ds)
    
    # subsurface mass
    #
    # maybe use a threshold of say -0.1 m?
    ds_sub=ds.where(ds.z<0.0) 
    # remove stranded data (in case subsurface particles intersected with land - it can happen especially with hz diffusion)
    ds_sub=ds_sub.where(~(ds_sub.status==stranded_flag))
    sub_mass=ds_sub.mass_oil.sum(dim='trajectory').rename('subsurface')

    # stranded mass
    #
    ds_strand=ds.where(ds.status==stranded_flag)
    strand_mass=ds_strand.mass_oil.sum(dim='trajectory').rename('stranded')
    
    # surface mass
    #
    # maybe use a threshold of say -0.1 m?
    ds_surf=ds.where(ds.z==0.0) # could use drop=True to reduce the size of the data, but should make no difference to end result
    # remove stranded data
    ds_surf=ds_surf.where(~(ds_surf.status==stranded_flag)) 
    surf_mass=ds_surf.mass_oil.sum(dim='trajectory').rename('surface')
    
    evap_mass=ds.mass_evaporated.sum(dim='trajectory').rename('evaporated')
    
    budget=xr.merge([evap_mass, surf_mass, sub_mass, strand_mass])
    budget.attrs["units"] = 'kg'
    
    budget.to_netcdf('trajectories_oil_budget.nc')
    
# if __name__ == "__main__":
    
        
