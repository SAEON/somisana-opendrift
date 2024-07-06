#!/usr/bin/env python
"""
library for computing the summary statistics from a set of stochastic simulations
"""

import sys, os, glob
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import opendrift
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.openoil import OpenOil
from opendrift.readers import reader_global_landmask
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from copy import copy
import xarray as xr
from xhistogram.xarray import histogram
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from bottleneck import push

class stochastic_mass_balance:
    """object used to compute the stochastic mass balance from stochastic iterations
    """
    def __init__(self, out_dir,
                 run_dir, 
                 date_start, 
                 run_id, 
                 run_id_end, 
                 increment_days, 
                 max_origin_marker=0
    ):
        self.out_dir                = out_dir # dir where output will be written
        self.run_dir                = run_dir # root dir where stochastic iterations were run
        self.date_start             = date_start
        self.run_id                 = run_id
        self.run_id_end             = run_id_end
        self.increment_days         = increment_days
        self.num_it                 = 0 # number of iterations initialised as zero
        self.max_origin_marker      = max_origin_marker
        
        # use a template file to initialise the 2D count where threshold is exceeded
        date_now = date_start+timedelta(days=(run_id-1)*increment_days)
        template_dir = run_dir+'run'+str(run_id).zfill(3)+'_'+date_now.strftime("%Y%m%d_%H%M%S")+'/'
        template_file = template_dir+'oil_budget_om'+str(max_origin_marker)+'.nc'
        ds = xr.open_dataset(template_file)
        days_since_start = (ds.time.data-ds.time.data[0]).astype('timedelta64[s]').astype(np.int32)/3600/24 # seriously convoluted but it works. Timedelta64 is not ideals
        ds = ds.assign(time = days_since_start)
        self.budget = ds
        
    def update(self, run_id):
        self.run_id = run_id
        date_now = self.date_start+timedelta(days=(run_id-1)*self.increment_days)
        # update the number of iterations
        self.num_it += 1 # count how many iterations have gone into computing the stochastic mass balance
        # add data for generating the stochastic budget
        if self.num_it>1:
            # but only if we're not on the first iteration (that was already added in __ini__)
            iteration_dir = self.run_dir+'run'+str(run_id).zfill(3)+'_'+date_now.strftime("%Y%m%d_%H%M%S")+'/'
            iteration_file = iteration_dir+'oil_budget_om'+str(self.max_origin_marker)+'.nc'
            ds = xr.open_dataset(iteration_file)
            days_since_start = (ds.time.data-ds.time.data[0]).astype('timedelta64[s]').astype(np.int32)/3600/24 # seriously convoluted but it works. Timedelta64 is not ideals
            ds = ds.assign(time = days_since_start)
            self.budget=xr.concat([self.budget, ds], dim="iteration")

    def __to_netcdf__(self):
        fname = self.out_dir+'oil_budget_om'+str(self.max_origin_marker)+'.nc'
        self.budget.to_netcdf(fname)

class stochastic_prob:
    """object used to compute the probability of oiling over defined threshold from stochastic iterations
    """
    def __init__(self, out_dir, 
                 run_dir, 
                 date_start, 
                 run_id, 
                 run_id_end, 
                 increment_days, 
                 filename_gridded,
                 threshold = 1.,
    ):
        self.out_dir                = out_dir # dir where output will be written
        self.run_dir                = run_dir # root dir where stochastic iterations were run
        self.date_start             = date_start
        self.run_id                 = run_id
        self.run_id_end             = run_id_end
        self.increment_days         = increment_days
        self.threshold              = threshold # units depends on posttype
        self.num_it                 = 0 # number of iterations initialised as zero
        self.posttype               = posttype
        self.filename_gridded       = filename_gridded 
        
        # use a template file to initialise the 2D count where threshold is exceeded
        date_now = date_start+timedelta(days=(run_id-1)*increment_days)
        template_dir = run_dir+'run'+str(run_id).zfill(3)+'_'+date_now.strftime("%Y%m%d_%H%M%S")+'/'
        template_file = template_dir+self.filename_gridded
        ds_template = xr.open_dataset(template_file)
        lon_bin = ds_template.lon_bin.data
        lat_bin = ds_template.lat_bin.data
        self.lon_bin = lon_bin
        self.lat_bin = lat_bin
        
        # initialise the dataset used for counting exceedance of threshold
        data_zeros=np.zeros((len(lon_bin), len(lat_bin)))
        data_nans=data_zeros.copy(); data_nans[data_nans==0]=np.nan
        self.stats = xr.Dataset(
            data_vars=dict(
                count_exceed=(["lat_bin","lon_bin"], data_zeros, {'units':'Number_of_times_threshold_is_exceeded'}),
                probability=(["lat_bin","lon_bin"], data_zeros, {'units':'Probability_threshold_is_exceeded'}),
                maximum=(["lat_bin","lon_bin"], data_zeros, {'units':'Maximum_over_all_runs'}),
                run_id_of_max=(["lat_bin","lon_bin"], data_nans, {'units':'Run_id_of_maximum'}),
                minimum_time=(["lat_bin","lon_bin"], data_zeros+10e6, {'units':'Minimum_time_over_all_runs'}),
            ),
            coords=dict(
                lon_bin=(["lon_bin"],lon_bin),
                lat_bin=(["lat_bin"],lat_bin)
            )
        )
        self.no_grid_cells_over_threshold = []
        self.stats.attrs["no_grid_cells_over_threshold"] = []
        self.run_ids = []
        self.stats.attrs["run_ids"] = []
        
    def update_prob(self, run_id):
        '''
        read the output from one of the deterministic runs and update the probability and minimum time variables
        '''
        self.run_id = run_id
        date_now = self.date_start+timedelta(days=(run_id-1)*self.increment_days)
        # update the number of iterations
        self.num_it += 1 # count how many iterations have gone into computing the probability
        # read data from this file
        iteration_dir = self.run_dir+'run'+str(run_id).zfill(3)+'_'+date_now.strftime("%Y%m%d_%H%M%S")+'/'
        iteration_file = iteration_dir+self.filename_gridded
        ds = xr.open_dataset(iteration_file)
        iteration_max = ds.maximum.data # 'maximum' is the standard variable name for all posttypes
        # update maximum if data from this iteration is greater
        previous_max=self.stats.maximum.data
        self.stats.maximum.data=np.maximum(iteration_max,previous_max)
        # save the run_id where this iteration is the overall maximum
        run_id_max=self.stats.run_id_of_max.data
        run_id_max[iteration_max>previous_max]=run_id
        self.stats.run_id_of_max.data=run_id_max
        # check where iteration max exceeds threshold and assign a value of 1
        iteration_exceed = np.zeros_like(iteration_max)
        iteration_exceed[iteration_max>self.threshold]=1
        # add up all the grid cells where the threshold is exceeded
        if len(self.no_grid_cells_over_threshold)==0:
            self.no_grid_cells_over_threshold=np.atleast_1d(np.sum(iteration_exceed))
            self.run_ids=np.atleast_1d(run_id)
        else:
            self.no_grid_cells_over_threshold=np.concatenate((self.no_grid_cells_over_threshold,np.atleast_1d(np.sum(iteration_exceed))),axis=0)
            self.run_ids=np.concatenate((self.run_ids,np.atleast_1d(run_id)),axis=0)
        self.stats.attrs["no_grid_cells_over_threshold"] = self.no_grid_cells_over_threshold # so that it is accessible in the output netcdf file we save later
        self.stats.attrs["run_ids"] =  self.run_ids
        # update count where threshold is exceeded
        self.stats.count_exceed.data += iteration_exceed
        # update probability of exceedance
        prob=self.stats.count_exceed.data/self.num_it
        # prob[prob==0]=np.nan
        self.stats.probability.data = prob
        
        # update the minimum time
        iteration_minimum_time = xr.where(iteration_exceed==1,ds.minimum_time.data,10e6) 
        # update minimum time if data from this iteration is smaller
        previous_minimum_time=self.stats.minimum_time.data
        self.stats.minimum_time.data=np.minimum(iteration_minimum_time,previous_minimum_time)
    
    def __to_netcdf__(self):
        fname = self.out_dir+self.filename_gridded.split('.')[0]+'_threshold'+str(self.threshold)
        self.stats.to_netcdf(fname + '.nc')
    
#if __name__ == "__main__":
    
    
