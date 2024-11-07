"""
some functions for processing raw opendrift output

"""

import os
import numpy as np
import xarray as xr
import opendrift
import xarray as xr
from xhistogram.xarray import histogram
import glob

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

def extents_from_lonlat(lon,lat,factor=0.01):
    '''
    factor = factor of domain size used to get dl, the buffer around the provided lon,lat
    '''
    lon_min = min(np.ravel(lon))
    lon_max = max(np.ravel(lon))
    lat_min = min(np.ravel(lat))
    lat_max = max(np.ravel(lat))
    dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
    extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    return extents

def grid_particles(fname,fname_out,grid_type='density',extents=None,dx_m=None,max_only=False):
    '''
    compute a eulerian particle density map from the output of an opendrift simulation
    fname = the file to do the gridding on
    fname_out = output netcdf filename
    grid_type = what kind of gridding to do. Options are 'density', 'surface_oil', 'stranded_oil'
    extents = the spatial extent of the grid [lon0,lon1,lat0,lat1]. If None, then this is automatically determined from the geographic extent of the particles
    dx_m = grid size in meters, if None, then a 50 x 50 regular grid is generated
    max_only = option to only write the maximum over the entire file to save disk space (boolean)
    '''
    
    def get_time_min(h):
        '''
        sub-function to extract the minimum time to occurence of a substance/object, 
        '''
        # compute the minimum time to oiling
        days_since_start = (h.time.data-h.time.data[0]).astype('timedelta64[s]').astype(np.int32)/3600/24 # seriously convoluted but it works. Timedelta64 is not ideals
        h_time=xr.full_like(h, 10e6) # using 10e6 as an arbitrary large number
        # loop through time and replace valid grid cells with time in days since spill start
        for ii, days in enumerate(days_since_start):
            h_time.data[ii,:,:]=xr.where(h.data[ii,:,:]>0,days,10e6)
        
        # compute the minimum time
        h_time_min=h_time.min(('time'))
        h_time_min = h_time_min.where(h_time_min != 10e6)
        h_time_min=h_time_min.rename('minimum_time')
        h_time_min.attrs["standard_name"] = 'minimum_time_to_occurrence'
        h_time_min.attrs["units"] = '_days_' # using underscores to avoid xarray automatically reading this variable as a timedelta64 later
        
        return h_time_min
    
    # get the data
    ds = xr.open_dataset(fname)
    ds = fill_deactivated(ds)
    lon = ds.lon.values
    lat = ds.lat.values
    
    if extents is None:
        extents = extents_from_lonlat(lon,lat)
    
    if dx_m is None:
        num_x = 50 # default number of grid points in each direction
        lonbin = np.linspace(extents[0], extents[1], num=num_x)
        latbin = np.linspace(extents[2], extents[3], num=num_x)
    else:
        lonbin, latbin = get_lonlat_bins(extents,dx_m)
    
    if grid_type == 'density':
        # compute the histogram of the particle locations
        # this will add up all the particles in each bin
        h = histogram(ds.lon,
                      ds.lat,
                      bins=[lonbin, latbin],
                      dim=['trajectory'],
                      density=False)
        # transpose to standard time,lat,lon ordering
        h = h.transpose("time", "lat_bin", "lon_bin")
        # normalise so that the data represents the particle density
        # i.e. the sum over the heat map should add up to 1
        # then we multiply by 100 to get percent
        h = 100 * h / len(lon)
        
        # convert oil volume per grid cell into oil thickness in micron
        h=h.rename('particle_density')
        h.attrs["units"] = '% particles'
        
        h_time_min = get_time_min(h)
        
        # compute the maximum density over the run
        h_max=h.max(('time'))
        h_max=h_max.rename('maximum')
        h_max.attrs["standard_name"] = 'maximum_particle_density'
        h_max.attrs["units"] = '% particles'
        
        if max_only:
            ds_out = xr.Dataset({'minimum_time': h_time_min, 'maximum': h_max})
            ds_out.to_netcdf(fname_out)
        else:
            ds_out = xr.Dataset({'particle_density': h, 'minimum_time': h_time_min, 'maximum': h_max})
            ds_out.to_netcdf(fname_out)
        
    elif grid_type == 'surface_oil':
    
        # identify stranded particles
        stranded_flag = get_stranded_flag(ds)
    
        # get data for computing surface thickness
        ds_surf=ds.where(ds.z==0.0)
        
        # remove stranded data (grid_stranded_oil() handles the stranded oil separately)
        ds_surf=ds_surf.where(~(ds_surf.status==stranded_flag)) 
        
        # compute the oil volume (m3) per surface particle
        # note we are using mass of the oil emulsion into account (see openoil.get_oil_budget)
        # mass_emulsion = mass_oil / (1 - water_fraction)
        surf_volume=(ds_surf.mass_oil / (1 - ds_surf.water_fraction))/ds_surf.density  
        
        # compute the histogram of surface volume
        # histogram omputes the number of particles per grid cell
        # setting 'weights=surf_volume' converts this to volume of surface oil per grid cell
        h_surf = histogram(ds_surf.lon,
                      ds_surf.lat,
                      bins=[lonbin, latbin],
                      dim=['trajectory'],
                      weights=surf_volume,
                      density=False)
        # transpose to standard time,lat,lon ordering
        h_surf = h_surf.transpose("time", "lat_bin", "lon_bin")
        
        # convert oil volume per grid cell into oil thickness in micron
        h_surf = h_surf/(np.power(dx_m,2))*1e6
        h_surf.attrs["standard_name"] = 'surface_oil_thickness'
        h_surf.attrs["units"] = 'micron'
        
        h_surf_time_min = get_time_min(h_surf)
        
        # compute the maximum surface thickness over the run
        h_surf_max=h_surf.max(('time'))
        h_surf_max=h_surf_max.rename('maximum')
        h_surf_max.attrs["standard_name"] = 'maximum_surface_oil_thickness'
        h_surf_max.attrs["units"] = 'micron'
        
        if max_only:
            ds_out = xr.Dataset({'minimum_time': h_surf_time_min, 'maximum': h_surf_max})
            ds_out.to_netcdf(fname_out)
        else:
            ds_out = xr.Dataset({'surface_thickness': h_surf, 'minimum_time': h_surf_time_min, 'maximum': h_surf_max})
            ds_out.to_netcdf(fname_out)
    
    elif grid_type == 'stranded_oil':
    
        # identify stranded particles
        stranded_flag = get_stranded_flag(ds)
    
        # get data for computing stranded oil concentration
        ds_strand=ds.where(ds.status==stranded_flag)
        
        # get the stranded mass per particle
        strand_mass=(ds_strand.mass_oil / (1 - ds_strand.water_fraction)) 
        
        # compute the histogram of stranded oil
        # histogram omputes the number of particles per grid cell
        # setting 'weights=strand_mass' converts this to the mass of stranded oil per grid cell
        h_strand = histogram(ds_strand.lon,
                      ds_strand.lat,
                      bins=[lonbin, latbin],
                      dim=['trajectory'],
                      weights=strand_mass,
                      density=False)
        # transpose to standard time,lat,lon ordering
        h_strand = h_strand.transpose("time", "lat_bin", "lon_bin")
        
        # convert mass of oil per grid cell from kg into g/m2
        beach_width_m=30 # m 1.5 m tidal range and 1:20 beach slope
        h_strand = h_strand*1000/dx_m/beach_width_m
        h_strand.attrs["standard_name"] = 'stranded_oil_concentration'
        h_strand.attrs["units"] = 'g/m2'
        
        h_strand_time_min = get_time_min(h_strand)
        
        # compute the maximum stranded concentration over the run
        h_strand_max=h_strand.max(('time'))
        h_strand_max=h_strand_max.rename('maximum')
        h_strand_max.attrs["standard_name"] = 'maximum_stranded_oil_concentration'
        h_strand_max.attrs["units"] = 'g/m2'
        
        if max_only:
            ds_out = xr.Dataset({'minimum_time': h_strand_time_min,'maximum': h_strand_max})
            ds_out.to_netcdf(fname_out)
        else:
            ds_out = xr.Dataset({'stranded_oil':h_strand, 'minimum_time': h_strand_time_min,'maximum': h_strand_max})
            ds_out.to_netcdf(fname_out)
            
    else:
        print('grid_type input of '+grid_type+' is not supported')
        print('valid options are density, surface_oil and stranded_oil')


def combine_gridded(dir_with_dirs,dirs=None,fname_gridded='gridded.nc',var_name=None):
    '''
    combine gridded output from multiple simulations into a single file
    the average values are computed over a selection of input gridded files, apart from minimum_time, in which case the minimum value is computed

    dir_with_dirs - path to a directory containing multiple directories, where each directory contains opendrift output
    dirs - a list of dir names inside dir_with_dirs which contain opendrift ouput. If None, then dirs is found dynamically from all dirs inside dir_with_dirs (except for the dirname 'combined')
    fname_gridded - the name of the gridded file which will be located inside each dir in dirs
    var_name - the time-dependent variable in fname_gridded. If None, then only the time-independent 'maximum' and 'minimum_time' variables are used
    '''

    if dirs is None:
        # find a list of all dirs inside dir_with_dirs, except for any dir called 'combined'
        all_dirs = glob.glob(dir_with_dirs + "/*/")  # Get all directories
        dirs = [d for d in all_dirs if not d.endswith("combined/")]  # Filter out 'combined' directories
    
    for i,dir_i in enumerate(dirs):
        fname_i = os.path.join(dir_with_dirs,dir_i,fname_gridded)
        ds = xr.open_dataset(fname_i)
        # maximum and minimum_time variables have dimensions (lat_bin, lon_bin) i.e. time independent
        maximum_i=ds.maximum
        minimum_time_i=ds.minimum_time
        if var_name is not None:
            var_i=ds[var_name]
        if i == 0:
            # initialise the output dataarrays using the first file
            maximum_out=maximum_i
            minimum_time_out=minimum_time_i
            if var_name is not None:
                var_out=var_i
        else:
            minimum_time_out = np.fmin(minimum_time_out, minimum_time_i)
            maximum_out = maximum_i + maximum_out
            if var_name is not None:
                # the var_name variable will have dimensions (time, lat_bin, lon_bin) i.e. time independent
                # But the length of the time dimension may not be consistent between different files
                # So we have to firstly align the two dataarrays along the time dimension
                var_out, var_i = xr.align(var_out, var_i, join="outer")
                
                # Fill missing times with zeros
                var_out = var_out.fillna(0)
                var_i = var_i.fillna(0)
                
                # Add the two DataArrays
                var_out = var_out + var_i

    # compute the averages for the maximum and var_name variables
    maximum_out=maximum_out/len(dirs)
    if var_name is not None:
        var_out=var_out/len(dirs)

    # make the 'combined' output directory
    combined_dir = os.path.join(dir_with_dirs,'combined')
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)
    
    # write the output file
    fname_out = os.path.join(combined_dir,fname_gridded)
    if var_name is not None:
        ds_out = xr.Dataset({var_name: var_out, 'minimum_time': minimum_time_out,'maximum': maximum_out})
        ds_out.to_netcdf(fname_out)
    else:
        ds_out = xr.Dataset({'minimum_time': minimum_time_out,'maximum': maximum_out})
        ds_out.to_netcdf(fname_out)

def combine_trajectories(dir_with_dirs,dirs=None,fname_traj='trajectories.nc'):
    '''
    combine raw opendrift output from multiple simulations into a single file
    by adding a new dimension called dir_name
        
    dir_with_dirs - path to a directory containing multiple directories, where each directory contains opendrift output
    dirs - a list of dir names inside dir_with_dirs which contain opendrift ouput. If None, then dirs is found dynamically from all dirs inside dir_with_dirs (except for the dirname 'combined')
    fname_traj - the name of the trajectories file (raw output) which will be located inside each dir in dirs
    '''
    
    if dirs is None:
        # find a list of all dirs inside dir_with_dirs, except for any dir called 'combined'
        all_dirs = glob.glob(dir_with_dirs + "/*/")  # Get all directories
        dirs = [d for d in all_dirs if not d.endswith("combined/")]  # Filter out 'combined' directories
    
    # I'm commenting this for now since I think it makes more sense to append
    # the data to the 'trajectory' dimension rather than create a new 'dir_name' dimension
    # But I'm not sure if there was a good reason for doing it like this in the first place?
    # datasets = []
    # for i,dir_i in enumerate(dirs):
    #     fname_i = os.path.join(dir_with_dirs,dir_i,fname_traj)
    #     ds = xr.open_dataset(fname_i)
    #     # Add a new dimension 'filename' to the dataset, using the file name as the label
    #     ds = ds.expand_dims({'dir_name': [dir_i.split('/')[-2]]})
    #     datasets.append(ds)
    # Concatenate datasets along the new 'filename' dimension
    # combined_ds = xr.concat(datasets, dim='dir_name')
    
    datasets = []
    for i, dir_i in enumerate(dirs):
        fname_i = os.path.join(dir_with_dirs, dir_i, fname_traj)
        ds = xr.open_dataset(fname_i)
        datasets.append(ds)
    
    # Concatenate datasets along the 'trajectory' dimension
    # the join='outer' argument should handle variable size in the time dimension
    combined_ds = xr.concat(datasets, dim='trajectory', join='outer')
    
    # make the 'combined' output directory (if not already there)
    combined_dir = os.path.join(dir_with_dirs,'combined')
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)
    
    # write the output file
    fname_out = os.path.join(combined_dir,fname_traj)
    combined_ds.to_netcdf(fname_out)

def oil_massbal(fname,fname_out):
    
    # # read the data
    oa = opendrift.open_xarray(fname)
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
    
    # evaporated
    evap_mass=ds.mass_evaporated.sum(dim='trajectory').rename('evaporated')
    
    budget=xr.merge([evap_mass, surf_mass, sub_mass, strand_mass])
    budget.attrs["units"] = 'kg'
    
    budget.to_netcdf(fname_out)
    
# if __name__ == "__main__":
    
        
