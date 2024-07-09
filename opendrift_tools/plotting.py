"""
functions for plotting / animating opendrift output
(there are some plotting methods as part of opendrift source code
 but we'd like to do our own tailor made ones for more flexibility)
"""

import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import cartopy.crs as ccrs
import pandas as pd
import opendrift_tools.postprocess as post
import matplotlib.path as mplPath
from opendrift.readers import reader_global_landmask

def get_croco_boundary(fname):
    '''
    Return lon,lat of perimeter around a CROCO grid (i.e. coordinates of bounding cells)
    '''
    with xr.open_dataset(fname) as ds:
        lon_rho = ds.lon_rho.values
        lat_rho = ds.lat_rho.values

    lon = np.hstack((lon_rho[0:, 0], lon_rho[-1, 1:-1],
                     lon_rho[-1::-1, -1], lon_rho[0, -2::-1]))
    lat = np.hstack((lat_rho[0:, 0], lat_rho[-1, 1:-1],
                     lat_rho[-1::-1, -1], lat_rho[0, -2::-1]))
    return lon, lat

def extents_2_polygon(extents):
    extents_poly = mplPath.Path(np.array([[extents[0], extents[2]],
                                    [extents[1], extents[2]],
                                    [extents[1], extents[3]],
                                    [extents[0], extents[3]]]))
    return extents_poly

def setup_plot(ax, lon, lat, extents=[], lscale='auto'):
    '''
    generic stuff applicable to all 2D plots
    lon = array of longitudes being plotted
    lat = array of latitudes being plotted
    extents = [lon_min, lon_max, lat_min, lat_max]
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
    '''
    # first need to get the domain extents if it's not set autmatically
    if len(extents) == 0:
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 # factor of domain size used to get dl
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    else:
        lon_min=extents[0]
        lon_max=extents[1]
        lat_min=extents[2]
        lat_max=extents[3]
    ax.set_extent(extents)
    
    # using opendrifts landmask plotting routine...
    reader_global_landmask.plot_land(ax, lon_min, lat_min, lon_max,
                                                 lat_max, False, lscale = lscale)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='dimgrey', alpha=0.5, linestyle=':')
    gl.right_labels = False
    gl.top_labels = False

def get_time_txt(ax,
                 time_plot,
                 time_start
                 ):
    '''
    get the time text to add to the plot, including 
    how long this time-step is from the start of the run
    '''
    days_since_start = (time_plot-time_start).astype('timedelta64[s]').astype(np.int32)/3600/24 # days
    days_since_start_int = np.floor(days_since_start).astype('int')
    hours_since_start = np.round((days_since_start-days_since_start_int)*24,0).astype('int')
    # time_plot=time_plot+np.timedelta64(2,'h') # show local time UTC+2 in plot?
    time_str = str(time_plot)[:10]+' '+str(time_plot)[11:19] # hack to get rid of the 'T' in the output of str(time_plot) - np.datetime64 is not great but is what our output comes as
    tx_time = time_str+'\n'+'('+str(days_since_start_int)+' days, '+str(hours_since_start)+' hours after start of release)'
    
    return tx_time

def add_text(ax,
                 tx,
                 loc=[0.5,1.01],
                 ):
    '''
    add the time text to the plot
    '''
    time_plt = ax.text(loc[0], loc[1], tx,
             ha='center', fontsize=12,
             transform=ax.transAxes)
    
    return time_plt

def add_cbar(var_plt,
             ticks=[],
             tick_font = 12,
             label='values',
             label_font=14,
             loc=[1., 0.2, 0.02, 0.6], # [left, bottom, width, height]
             orientation='vertical'):
    
    '''
    Add a colorbar to a plot
    '''
    
    cbarax = plt.gcf().add_axes(loc) 
    cbar_plt = plt.colorbar(var_plt, cbarax,
                        ticks=ticks,
                        orientation=orientation)
    cbar_plt.set_label(label, fontsize=label_font)
    cbar_plt.ax.tick_params(labelsize=tick_font)
    
    return cbar_plt

def plot_particles(fname,
        var_str='z', # variable to plot
        tstep=0, # the step to plot, or the first step to animate.
        # options related to the figure layout
        figsize=(6,6), # (hz,vt)        
        extents = [], # [lon0,lon1,lat0,lat1]
        lscale = 'auto', # resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
        # options relating to the release location
        lon_release=None, 
        lat_release=None,
        size_release = 50,
        # options relating to the dispaly of data, colormap and colorbar
        size_scat = 20, # size of the scatter data to be plotted
        ticks = [-9999,-9998], # the ticks to plot relating to the colormap (can be irregularly spaced)
        cmap = 'Greys', # colormap to use
        plot_cbar = False,
        cbar_loc = [0.9, 0.2, 0.02, 0.6], # where on the plot to put the colorbar
        cbar_label = 'depth (m)',
        # options related to the plot output file
        jpg_out=None, # filename of the jpg file
        write_jpg=False,
        # options related to the animation
        gif_out=None, # filename of the gif file
        write_gif=False,
        skip_time = 1, # every nth time-step will be animated (if provided)
        tstep_end=None, # The last timestep to animate. Only used if write_gif = True.
        ):
    '''
    this is a convenience function for doing a quick 2D plot of particles with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    # get the data
    ds = xr.open_dataset(fname)
    ds = post.fill_deactivated(ds)
    time_start = ds.time.data[0]
    
    # subset to the time-step
    ds_tstep = post.subset_tstep(ds,tstep)
    time_plot = ds_tstep.time.values
    lon=ds_tstep.lon.values
    lat=ds_tstep.lat.values
    var_data = ds_tstep[var_str].values
    
    # identify stranded particles so we can plot them as red
    stranded_flag = post.get_stranded_flag(ds)
    ds_strand_tstep=ds_tstep.where(ds_tstep.status==stranded_flag)
    lon_strand=ds_strand_tstep.lon.values
    lat_strand=ds_strand_tstep.lat.values
    
    # set up the plot
    fig = plt.figure(figsize=figsize) 
    ax = plt.axes(projection=ccrs.Mercator())
    setup_plot(ax,lon,lat,extents=extents,lscale=lscale)
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    # add the particles
    scat = ax.scatter(lon,lat, s=size_scat, 
                      c=var_data,
                      cmap=cmap,
                      norm=cmap_norm,
                      transform=ccrs.PlateCarree())
    
    # add the stranded particles (hard coded as red)
    scat_strand = ax.scatter(lon_strand,lat_strand, s=size_scat, color='r', 
                      transform=ccrs.PlateCarree())
    
    if lon_release is not None:
        ax.scatter(lon_release,lat_release, size_release, transform=ccrs.PlateCarree(),marker='X',color='k')
    
    tx_time = get_time_txt(ax, time_plot, time_start)
    time_plt = add_text(ax,tx_time,loc=[0.5,1.01])
    
    if plot_cbar:
        add_cbar(scat,label=cbar_label,ticks=ticks,loc=cbar_loc)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            # automatically come up with a file name
            jpg_out = fname.split('.nc')[0]+'_'+var_str+'_'+pd.to_datetime(time_plot).strftime("%Y%m%d_%H")+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    # and/or write a gif if specified
    if write_gif: # do the animation
        def plot_tstep(i):
            
            ds_tstep = post.subset_tstep(ds,i)
            time_plot = ds_tstep.time.values
            lon=ds_tstep.lon.values
            lat=ds_tstep.lat.values
            var_data = ds_tstep[var_str].values
            
            ds_strand_tstep=ds_tstep.where(ds_tstep.status==stranded_flag)
            lon_strand=ds_strand_tstep.lon.values
            lat_strand=ds_strand_tstep.lat.values
            
            # update the time label
            tx_time = get_time_txt(ax, time_plot, time_start)
            time_plt.set_text(tx_time)
            
            # update the scatter data
            scat.set_offsets(np.c_[lon,lat])
            scat.set_array(var_data)
            scat_strand.set_offsets(np.c_[lon_strand,lat_strand])
        
        # animate
        if tstep_end is None: # if not defined then animate to end of file 
            tstep_end = len(ds.time) - 1
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(tstep,tstep_end,skip_time)) 
        
        if gif_out is None:
            gif_out = fname.split('.nc')[0]+'_'+var_str+'.gif'
        anim.save(gif_out, writer='imagemagick')

def lonlat_2_corners(lon,lat):
    '''
    take lon lat data from gridded output and return new lon,lat data which
    represent the grid cell corners instead of the centre of the grid cell
    '''
    # Calculate the spacing between original lon and lat points
    lon_diff = np.diff(lon)
    lat_diff = np.diff(lat)
      
    # Create empty grids with one extra element in each dimension
    lon_out = np.zeros((len(lon)+1,), dtype=float)
    lat_out = np.zeros((len(lat)+1,), dtype=float)
    
    # Fill the interior grid points with averages of original coordinates
    lon_out[1:-1] = lon[:-1] + lon_diff/2
    lat_out[1:-1] = lat[:-1] + lat_diff/2
    
    # add the end points
    lon_out[0] = lon[0] - lon_diff[0]/2
    lat_out[0] = lat[0] - lat_diff[0]/2
    lon_out[-1] = lon[-1] + lon_diff[-1]/2
    lat_out[-1] = lat[-1] + lat_diff[-1]/2
    
    return lon_out,lat_out

def plot_gridded(fname,
        var_str='particle_density', # variable to plot
        tstep=0, # the step to plot, or the first step to animate.
        # options related to the figure layout
        figsize=(6,6), # (hz,vt)        
        extents = [], # [lon0,lon1,lat0,lat1]
        lscale = 'auto', # resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
        # options relating to the release location
        lon_release=None, 
        lat_release=None,
        size_release = 50,
        # options relating to the dispaly of data, colormap and colorbar
        ticks = np.linspace(0,0.2,num=11), # the ticks to plot relating to the colormap (can be irregularly spaced)
        cmap = 'Spectral_r', # colormap to use
        plot_cbar = False,
        cbar_loc = [0.9, 0.2, 0.02, 0.6], # where on the plot to put the colorbar
        cbar_label = 'particle density (-)',
        # options related to the plot output file
        jpg_out=None, # filename of the jpg file
        write_jpg=False,
        # options related to the animation
        gif_out=None, # filename of the gif file
        write_gif=False,
        skip_time = 1, # every nth time-step will be animated (if provided)
        tstep_end=None, # The last timestep to animate. Only used if write_gif = True.
        ):
    '''
    this is a convenience function for doing a quick 2D plot of gridded output with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    # get the data
    ds = xr.open_dataset(fname)
    time_start = ds.time.data[0]
    
    # only plot where data shows up
    ds=ds.where(ds[var_str]>0.0)
    
    # subset to the time-step
    ds_tstep = post.subset_tstep(ds,tstep)
    time_plot = ds_tstep.time.values
    lon=ds_tstep.lon_bin.values
    lat=ds_tstep.lat_bin.values
    lon,lat=lonlat_2_corners(lon,lat)
    var_data = ds_tstep[var_str].values
    
    # set up the plot
    fig = plt.figure(figsize=figsize) 
    ax = plt.axes(projection=ccrs.Mercator())
    setup_plot(ax,lon,lat,extents=extents,lscale=lscale)
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    # plot the data
    var_plt = ax.pcolormesh(lon,
                              lat,
                              var_data,
                              cmap=cmap,
                              norm=cmap_norm,
                              transform=ccrs.PlateCarree())
    
    if lon_release is not None:
        ax.scatter(lon_release,lat_release, size_release, transform=ccrs.PlateCarree(),marker='X',color='k')
    
    tx_time = get_time_txt(ax, time_plot, time_start)
    time_plt = add_text(ax,tx_time,loc=[0.5,1.01])
    
    if plot_cbar:
        add_cbar(var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            # automatically come up with a file name
            jpg_out = fname.split('.nc')[0]+'_'+var_str+'_'+pd.to_datetime(time_plot).strftime("%Y%m%d_%H")+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
    # and/or write a gif if specified
    if write_gif: # do the animation
        def plot_tstep(i):
            
            ds_tstep = post.subset_tstep(ds,i)
            time_plot = ds_tstep.time.values
            var_i = ds_tstep[var_str].values
            var_plt.set_array(var_i.ravel())
            
            # update the time label
            tx_time = get_time_txt(ax, time_plot, time_start)
            time_plt.set_text(tx_time)
            
        
        # animate
        if tstep_end is None: # if not defined then animate to end of file 
            tstep_end = len(ds.time) - 1
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(tstep,tstep_end,skip_time)) 
        
        if gif_out is None:
            gif_out = fname.split('.nc')[0]+'_'+var_str+'.gif'
        anim.save(gif_out, writer='imagemagick')

def plot_gridded_stats(fname,
        var_str='probability', # variable to plot
        # options related to the figure layout
        figsize=(6,6), # (hz,vt)        
        extents = [], # [lon0,lon1,lat0,lat1]
        lscale = 'auto', # resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
        # options relating to the release location
        lon_release=None, 
        lat_release=None,
        size_release = 50,
        # options relating to the dispaly of data, colormap and colorbar
        ticks = np.linspace(0,1,num=11), # the ticks to plot relating to the colormap (can be irregularly spaced)
        cmap = 'Spectral_r', # colormap to use
        plot_cbar = False,
        cbar_loc = [0.9, 0.2, 0.02, 0.6], # where on the plot to put the colorbar
        cbar_label = 'probability of occurrence (-)',
        # options related to the plot output file
        jpg_out=None, # filename of the jpg file
        write_jpg=False
        ):
    '''
    this is a convenience function for doing a quick 2D plot of summary stats gridded output with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    # get the data
    ds = xr.open_dataset(fname)
    
    # only plot where data shows up
    ds=ds.where(ds[var_str]>0.0)
    
    # subset to the time-step
    lon=ds.lon_bin.values
    lat=ds.lat_bin.values
    lon,lat=lonlat_2_corners(lon,lat)
    var_data = ds[var_str].values
    
    # set up the plot
    fig = plt.figure(figsize=figsize) 
    ax = plt.axes(projection=ccrs.Mercator())
    setup_plot(ax,lon,lat,extents=extents,lscale=lscale)
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    # plot the data
    var_plt = ax.pcolormesh(lon,
                              lat,
                              var_data,
                              cmap=cmap,
                              norm=cmap_norm,
                              transform=ccrs.PlateCarree())
    
    if lon_release is not None:
        ax.scatter(lon_release,lat_release, size_release, transform=ccrs.PlateCarree(),marker='X',color='k')
        
    if plot_cbar:
        add_cbar(var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc)
    
    # write a jpg if specified
    if write_jpg:
        if jpg_out is None:
            # automatically come up with a file name
            jpg_out = fname.split('.nc')[0]+'_'+var_str+'.jpg'
        plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
    
# if __name__ == "__main__":
    

    
