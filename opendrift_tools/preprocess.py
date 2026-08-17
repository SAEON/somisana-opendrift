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
import opendrift_tools.postprocess as post

def array_to_txt(array, filename):
    with open(filename, 'w') as file:
        for item in array:
            file.write(str(item) + '\n')

def array_from_txt(filename, dtype=float):
    """
    Read values from a text file into a NumPy array.

    Parameters
    ----------
    filename : str
        Path to text file.
    dtype : type
        Data type to convert values to.
        Examples: float, int, str

    Returns
    -------
    np.ndarray
    """

    with open(filename, 'r') as file:
        array = [dtype(line.strip()) for line in file]

    return np.array(array, dtype=dtype)

def get_particle_positions(flag_file,chl_file,config_dir,domain=None,Max_particles=5000):
    """
    Function to define particle positions within a bloom area based of phytoplankton flags
    and chlorophyll concentrations dessiminated by CSIR. 
    The distribution of the particles is largely dependant on the size of the 
    domain and the max number of particles. 
    
    Inputs:
        flag_file: CSIR phytoplankton flags file. We seed particles at grid cells
                  flagged 6 (i.e. Red tide).
        chl_file: CSIR Chlorophyll-a file. Max number of grid cells are distributed 
                between the grid cells flagged 6. The number of particles seeded within 
                the flags depends on the chl concentration.
        config_dir: Directory path of the config.py file used to set up the opendrift
                    experiment.
        domain: List of minimum and maximum longitudes and latitudes to subset 
                spatially [lon_min, lon_max, lat_min, lat_max].  
        Max_particles: Maximum number of particles to seed within the domain (int).
    
    Outputs:
        lon_rand: longitudes of the particles in a 1D array.
        lat_rand:  latitudes of the particles in a 1D array.
        In addition, we save the chlorophyll-a concentration associated to each particle 
        in a seperate .txt file which can be used to convert back to a Eularian field. 
    """
    try:
        # load the chl data and spatial dimensions from CSIR
        if domain is None:
            ds_chl = xr.open_dataset(chl_file).isel(band=0)
            ds_flags = xr.open_dataset(flag_file).isel(band=0)
        else:
            ds_chl  = xr.open_dataset(chl_file).isel(band=0
                                                     ).sel(x=slice(domain[0],domain[1]),
                                                           y=slice(domain[3],domain[2])
                                                           )
            ds_flags  = xr.open_dataset(flag_file).isel(band=0
                                                        ).sel(x=slice(domain[0],domain[1]),
                                                              y=slice(domain[3],domain[2])
                                                              )
        
        # Read chlorophyll data from CSIR
        chl_concentration = ds_chl.chl.values
        lon = ds_chl.x.values
        lat = ds_chl.y.values
        ds_chl.close()
        
        # read phytoplankton flag data from CSIR
        flags = ds_flags.phytoplankton.values
        LON,LAT = np.meshgrid(lon,lat)
        ds_flags.close()
        
        array_to_txt(lon,config_dir + '/lons.txt'), array_to_txt(lat,config_dir + '/lats.txt')

        # compute grid area
        grid_cell_area = post.compute_grid_area(lon, lat)

        # Mask all chlorophyll values that are not flagged 6
        chl_concentration[np.isnan(chl_concentration)] = 0
        chl_concentration[flags!=3] = 0
    
        # Compute the chlorophyll mass each grid cell
        chl_mass = chl_concentration * grid_cell_area
                
        # Compute the chlorophyll budget
        chl_budget = np.sum(chl_mass)
    
        # Calculate chlorophyll percentage per grid cell
        chl_percent = chl_mass / chl_budget  * 100
        
        # Calculate number of particles within each grid cell based on the max number of 
        # particles and the chlorophyll percentage of each grid cell
        Nparticles_gc = chl_percent / 100 * Max_particles 
        
        # We do round it since the number of particles within a grid cell cannot hava decimals.
        Nparticles_gc = np.asanyarray(np.floor(Nparticles_gc), dtype = int )
        
        # We get the new total number of particles that we will seed. 
        # This will be different to the original amount due to rounding.
        Nparticles = int(np.sum(Nparticles_gc))

        # Compute the chlorophyll mass associated with each particle. 
        # This number is important for when we regrid back to a Eularian field. 
        array_to_txt(np.array([chl_budget / Nparticles]),
                     config_dir + '/chl_mass.txt')
        
        # we get the indices of the gridcells which contains particles
        idx = np.argwhere(Nparticles_gc >= 1)
        j,i = idx[:,0],idx[:,1]
        
        # Make empty arrays to write the particle positions to
        lon_rand, lat_rand = np.zeros(Nparticles), np.zeros(Nparticles)
        i_start = 0
        for ii in range(i.size):
            n_points = Nparticles_gc[j[ii],i[ii]]
            i_end = i_start+n_points
            lon_rand[i_start:i_end], lat_rand[i_start:i_end] = random_points_in_cell(LON, 
                                                                                     LAT, 
                                                                                     i[ii], 
                                                                                     j[ii], 
                                                                                     n_points
                                                                                     )
            i_start=i_start+n_points

        return lon_rand,lat_rand
    
    except:
        return None,None
        

def random_points_in_cell(lon, lat, i, j, n_points):
    """
    Generate random (lon, lat) points uniformly inside a grid cell for
    cell-centered 2D grids (e.g. ROMS/CROCO lon_rho, lat_rho).

    Parameters
    ----------
    lon, lat : 2D arrays
        Longitude and latitude arrays at cell centers (same shape).
    i, j : int
        Indices of the target grid cell center.
        Must not be on the boundary (i > 0, j > 0, i < lon.shape[0]-1, j < lon.shape[1]-1)
    n_points : int
        Number of random points to generate inside the cell.

    Returns
    -------
    lon_rand, lat_rand : 1D arrays
        Random longitude and latitude values inside the cell.
    """

    # --- Compute approximate cell corners by averaging neighboring centers ---
    lon_nw = 0.25 * (lon[j, i] + lon[j-1, i] + lon[j, i-1] + lon[j-1, i-1])
    lon_ne = 0.25 * (lon[j, i] + lon[j-1, i] + lon[j, i+1] + lon[j-1, i+1])
    lon_se = 0.25 * (lon[j, i] + lon[j+1, i] + lon[j, i+1] + lon[j+1, i+1])
    lon_sw = 0.25 * (lon[j, i] + lon[j+1, i] + lon[j, i-1] + lon[j+1, i-1])

    lat_nw = 0.25 * (lat[j, i] + lat[j-1, i] + lat[j, i-1] + lat[j-1, i-1])
    lat_ne = 0.25 * (lat[j, i] + lat[j-1, i] + lat[j, i+1] + lat[j-1, i+1])
    lat_se = 0.25 * (lat[j, i] + lat[j+1, i] + lat[j, i+1] + lat[j+1, i+1])
    lat_sw = 0.25 * (lat[j, i] + lat[j+1, i] + lat[j, i-1] + lat[j+1, i-1])

    # --- Define corners of the grid cells ---
    lon_corners = [lon_nw, lon_ne, lon_se, lon_sw]
    lat_corners = [lat_nw, lat_ne, lat_se, lat_sw]
    
    # --- Helper: random points inside one triangle using barycentric coordinates ---
    def rand_in_triangle(p1, p2, p3, n):
        # P = ( 1 − r1 <200b>) * p1 <200b>+ r1 * <200b>( 1 − r2 <200b>) * p2 <200b>+ r1<200b> * r2 * <200b>p3<200b>
        # Take note that we square a random point (in this case r1) to ensure that we get 
        # a random uniform distribution throughout our triangle. 
        r1 = np.sqrt(np.random.rand(n))
        r2 = np.random.rand(n)
        lon_r = (1 - r1) * p1[0] + r1 * (1 - r2) * p2[0] + r1 * r2 * p3[0]
        lat_r = (1 - r1) * p1[1] + r1 * (1 - r2) * p2[1] + r1 * r2 * p3[1]
        return lon_r, lat_r

    # Split total points between the two triangles
    n1 = n_points // 2
    n2 = n_points - n1

    p0 = (lon_corners[0], lat_corners[0])  # NW
    p1 = (lon_corners[1], lat_corners[1])  # NE
    p2 = (lon_corners[2], lat_corners[2])  # SE
    p3 = (lon_corners[3], lat_corners[3])  # SW

    # Triangle 1: p0, p1, p2
    lon1, lat1 = rand_in_triangle(p0, p1, p2, n1)
    # Triangle 2: p0, p2, p3
    lon2, lat2 = rand_in_triangle(p0, p2, p3, n2)

    # Combine both sets of points
    lon_rand = np.concatenate([lon1, lon2])
    lat_rand = np.concatenate([lat1, lat2])

    return lon_rand, lat_rand

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
        
        if getattr(config, "use_croco_grid", False):
            croco_grids = config.croco_grid
            
            if len(croco_files) != len(croco_grids):
                raise ValueError("The number of CROCO data files and grid files must match.")

            # There seems to be some mismatch between variables
            # So we force the merge - probably not the best
            # We first store the original merge function
            original_merge = xr.merge 
            
            # Create a wrapper function that forces 'override'
            def override_merge(*args, **kwargs):
                kwargs['compat'] = 'override'
                return original_merge(*args, **kwargs)
                
            for croco_file, croco_grid in zip(croco_files, croco_grids):
                # Apply the patch so OpenDrift ignores the slight differences
                xr.merge = override_merge 
                
                reader_croco = reader_ROMS_native.Reader(croco_file, gridfile=croco_grid)
                
                # Revert the patch immediately so we don't break anything else!
                xr.merge = original_merge 
                
                croco_ref_time = datetime(config.croco_Yorig, 1, 1)
                reader_croco = set_croco_time(reader_croco, croco_ref_time)
                o.add_reader(reader_croco)
                print(f'\n We added CROCO successfully: {croco_file} \n')

                
            for croco_file, croco_grid in zip(croco_files, croco_grids):
                reader_croco = reader_ROMS_native.Reader(croco_file, gridfile=croco_grid)
                croco_ref_time = datetime(config.croco_Yorig, 1, 1)
                reader_croco = set_croco_time(reader_croco, croco_ref_time)
                o.add_reader(reader_croco)
                print(f'\n We added CROCO successfully: {croco_file} \n')
                
        else:
            # Standard loop if no grids are provided
            for croco_file in croco_files:
                reader_croco = reader_ROMS_native.Reader(croco_file)
                croco_ref_time = datetime(config.croco_Yorig, 1, 1)
                reader_croco = set_croco_time(reader_croco, croco_ref_time)
                o.add_reader(reader_croco)
                print(f'\n We added CROCO successfully: {croco_file} \n')
                
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
