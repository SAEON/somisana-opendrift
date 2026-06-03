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
import os

def write_seed_particles(
        outfile,
        lonp,
        latp,
        weightp,
        grid_lon=None,
        grid_lat=None,
        release_time=None):
    """
    Write seeded particles to NetCDF.

    Parameters
    ----------
    outfile : str
        Output NetCDF filename.

    lonp : ndarray
        Particle longitude.

    latp : ndarray
        Particle latitude.

    weightp : ndarray
        Particle weight (mg).
        
    grid_lon : ndarray (1D)
        longitudes of the grid

    grid_lat : ndarray (1D)
        latitudes of the grid
        
    release_time : datetime-like, optional
        Particle release time.
    """
    
    if os.path.exists(outfile):
        os.remove(outfile)

    n_particles = len(lonp)

    ds = xr.Dataset(
        data_vars={
            "lon": (
                ["particle"],
                np.asarray(lonp),
                {
                    "long_name": "particle longitude",
                    "units": "degrees_east"
                }
            ),

            "lat": (
                ["particle"],
                np.asarray(latp),
                {
                    "long_name": "particle latitude",
                    "units": "degrees_north"
                }
            ),

            "weight": (
                ["particle"],
                np.asarray(weightp),
                {
                    "long_name": "chlorophyll mass represented by particle",
                    "units": "mg"
                }
            )
        },

        coords={
            "particle": np.arange(n_particles)
        },

        attrs={
            "title": "Seeded phytoplankton particles",
            "Conventions": "CF-1.8"
        }
    )
    
    if grid_lon is not None:
        ds["grid_lon"] = xr.DataArray(
            np.array(grid_lon),
            dims=["x"],
            attrs={
                "long_name": "grid longitude",
                "units": "degrees_east"
            }
        )
            
    if grid_lat is not None:
        ds["grid_lat"] = xr.DataArray(
            np.array(grid_lat),
            dims=["y"],
            attrs={
                "long_name": "grid latitude",
                "units": "degrees_north"
            }
        )

    if release_time is not None:
        ds["release_time"] = xr.DataArray(
            np.repeat(
                np.datetime64(release_time),
                n_particles
            ),
            dims=["particle"]
        )
        
    ds.to_netcdf(outfile)

    print(f"Saved {n_particles} particles to:")
    print(outfile)

def seed_particles_from_chlorophyll(
        flag_file,
        chl_file,
        domain,
        dz=5,
        nmax=100,
        earth_radius=6371000.0):
    """
    Seed particles from a chlorophyll field using adaptive
    logarithmic particle densities.

    Parameters
    ----------
    flag_file : str
        Phytoplankton flag file.

    chl_file : str
        Chlorophyll concentration file.

    domain : list
        [xmin, xmax, ymax, ymin]

    dz : float
        Layer thickness (m).

    nmax : int
        Maximum particles per cell.

    Returns
    -------
    lons : ndarray
        Grid longitudes

    lats : ndarray
        Grid latitudes
        
    lonp : ndarray
        Particle longitudes.

    latp : ndarray
        Particle latitudes.

    weightp : ndarray
        Particle weights.

    Np : ndarray
        Number of particles per source cell.

    source_i : ndarray
        Source i-indices.

    source_j : ndarray
        Source j-indices.
    """

    # -------------------------------------------------
    # Get grid cells that are flagged 6 - i.e. Red tide
    # -------------------------------------------------

    j, i = get_flag_indexes(
        file_in=flag_file,
        domain=domain
    )

    if j is None or len(j) == 0:
        return (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([])
        )

    # -------------------------------------------------
    # Read chlorophyll concentration
    # -------------------------------------------------

    ds = xr.open_dataset(chl_file).sel(
        x=slice(domain[0], domain[1]),
        y=slice(domain[3], domain[2])
    )

    chl_conc = ds.chl.values.squeeze()

    lons = ds.x.values
    lats = ds.y.values

    LONS, LATS = np.meshgrid(lons, lats)

    ds.close()

    # -------------------------------------------------
    # Grid spacing
    # -------------------------------------------------

    dx, dy = compute_grid_spacing(
        LONS,
        LATS,
        R=earth_radius
    )

    # -------------------------------------------------
    # Particle numbers
    # -------------------------------------------------

    Np = compute_particle_number(
        chl_conc[j, i],nmax=nmax
    )

    # -------------------------------------------------
    # Vilidity mask
    # -------------------------------------------------

    C_seed = chl_conc[j, i]
    dx_seed = dx[j, i]
    dy_seed = dy[j, i]

    valid = (
        np.isfinite(C_seed) &
        np.isfinite(dx_seed) &
        np.isfinite(dy_seed) &
        np.isfinite(Np) &
        (Np > 0)
    )

    j = j[valid]
    i = i[valid]

    C_seed = C_seed[valid]
    dx_seed = dx_seed[valid]
    dy_seed = dy_seed[valid]

    Np = Np[valid]

    if len(Np) == 0:
        return (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([])
        )

    # -------------------------------------------------
    # Particle weights
    # -------------------------------------------------

    weights = compute_particle_weight(
        C_seed,
        dx_seed,
        dy_seed,
        dz,
        Np
    )

    # -------------------------------------------------
    # Create particles
    # -------------------------------------------------

    lonp = []
    latp = []
    weightp = []

    for n in range(len(Np)):

        lon_rand, lat_rand = random_points_in_cell(
            LONS,
            LATS,
            i[n],
            j[n],
            Np[n]
        )

        lonp.extend(lon_rand.tolist())
        latp.extend(lat_rand.tolist())

        weightp.extend(
            [weights[n]] * Np[n]
        )

    return lons,  lats, np.asarray(lonp), np.asarray(latp), np.asarray(weightp), Np, i, j

def get_flag_indexes(file_in,domain=None):
    #
    # Function to get the indexes at points that are flagged 
    # as a red tide in the phytoplankton flag file. 
    #
    try:
        if domain is not None:
            ds_flags = ds_flags = xr.open_dataset(file_in).sel(x=slice(domain[0],domain[1]),
                                                               y=slice(domain[3],domain[2])
                                                               )
        else:
            ds_flags = xr.open_dataset(file_in)
    except:
        print('\nError extracting flagged points. No flag data.\n')
        return None,None
        
    phytoplankton_flags = ds_flags.phytoplankton.values.squeeze()
    ds_flags.close()
    
    idx = np.argwhere((phytoplankton_flags==6))
    j,i = idx[:,0],idx[:,1]
    
    return j,i

def compute_particle_number(C,cmin=0.1,cmax=300,nmin=1,nmax=100):
    """
    Compute number of particles per grid cell to seed based on chlorophyll concentration. 
    A greater number of particles are assigned to a grid cell with a high 
    clorophyll concentration and fewer particles are assigned to a grid cell 
    with a lower clorophyll concentration.
    
    We use the following log function : 
        Np = nmin + (nmax - nmin) x (log(C / cmin) / log(cmax / cmin))
        
    Parameters
    ----------
    C : float or ndarray
        Chlorophyll concentration (mg m^-3)
    cmin : float 
        Minimum Chlorophyll concentration (mg m^-3)
    cmax : float 
        Maximum Chlorophyll concentration (mg m^-3)
    nmin : float 
        Minimum number of particles    
    nmax : float 
        Maximum number of particles    
        
    Returns
    -------
    Np : int or ndarray
        Number of particles to seed
    """
    
    # Small value to avoid log(0)
    eps=1e-12
    
    # Convert to array
    C = np.asarray(C, dtype=float)

    # Replace invalid values
    C = np.where(np.isfinite(C), C, eps)

    # Prevent zero or negative values
    C = np.maximum(C, eps)

    # Logarithmic scaling
    Np = (
        nmin
        + (nmax - nmin)
        * np.log(C / cmin)
        / np.log(cmax / cmin)
    )

    # Clip range
    Np = np.clip(Np, nmin, nmax)
    
    # round to nearest integer
    Np = np.rint(Np).astype(int)

    return Np

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

def compute_particle_weight(C,dx,dy,dz,Np):
    """
    Compute weight carried by each particle.

    Parameters
    ----------
    C : float or ndarray
        Chlorophyll concentration (mg m^-3)

    dx, dy, dz : float or ndarray
        Grid-cell dimensions (m)

    Np : int or ndarray
        Number of particles in cell

    Returns
    -------
    w : float or ndarray
        Weight per particle (mg)
    """

    # Cell volume
    V = dx * dy * dz

    # Total chlorophyll mass in cell
    M = C * V

    # Weight per particle
    w = M / Np

    return w

def compute_grid_spacing(lon, lat, R=6371000.0):
    """
    Compute dx and dy grid spacing from lon/lat cell centers.

    Parameters
    ----------
    lon : 2D ndarray
        Longitudes in degrees

    lat : 2D ndarray
        Latitudes in degrees

    Returns
    -------
    dx : 2D ndarray
        Zonal grid spacing (m)

    dy : 2D ndarray
        Meridional grid spacing (m)
    """
    
    # Convert to radians
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    
    dlon = np.gradient(lon_rad, axis=1)
    dlat = np.gradient(lat_rad, axis=0)

    # dx depends on latitude
    dx = R * np.cos(lat_rad) * dlon

    # dy
    dy = R * dlat

    return np.abs(dx), np.abs(dy)


def get_seed_points(file_in,domain=None):
    #
    # Function to get the longoitude and latitude positions at points that are flagged 
    # as a red tide in the phytoplankton flag file. 
    #
    try:
        if domain is not None:
            ds_flags = ds_flags = xr.open_dataset(file_in).sel(x=slice(domain[0],domain[1]),
                                                               y=slice(domain[3],domain[2])
                                                               )
        else:
            ds_flags = xr.open_dataset(file_in)
    except:
        print('\nError extracting flagged points. No flag data.\n')
        return None,None
        
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
