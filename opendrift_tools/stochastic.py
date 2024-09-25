#!/usr/bin/env python
"""
tools for carrying out stochastic simulations, doing some postprocessing and summary statistics
"""

import sys, os, shutil
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
from matplotlib.path import Path
from opendrift_tools.run import oil as run_oil
from opendrift_tools.run import leeway as run_leeway
from opendrift_tools.run import oceandrift as run_oceandrift
from opendrift_tools.postprocess import grid_particles, oil_massbal

class base_stochastic:
    '''
    base class for stochastic simulations and postprocessing - inherited by other derived classes
    '''
    def __init__(self,
                 run_dir, 
                 date_start, 
                 run_id,
                 increment_days,
                 run_id_end
    ):
        self.run_dir                = run_dir # root dir where stochastic iterations were run
        self.date_start             = date_start # datetime of first run
        self.run_id                 = run_id # run id to start on (you don't have to start at run001)
        self.increment_days         = increment_days # number of days increment between runs
        self.run_id_end             = run_id_end # run id to end on
        self.num_it                 = 1 # number of iterations initialised as one (used to compute statistics)
        self.update_iteration_dir()
        
    def update_iteration_dir(self):
        '''
        get the directory of the current stochastic iteration
        '''
        date_now = self.date_start+timedelta(days=(self.run_id-1)*self.increment_days)
        self.iteration_dir = self.run_dir+'/run'+str(self.run_id).zfill(3)+'_'+date_now.strftime("%Y%m%d_%H")+'/'
    
class run_stochastic(base_stochastic):
    '''
    class used to run a bunch of OpenDrift simulations as part of a stochastic analysis
    '''
    def __init__(self, run_dir, date_start, run_id, increment_days, run_id_end, model_type='oil'):
        super().__init__(run_dir, date_start, run_id, increment_days, run_id_end) # this is where the base_stochastic class is inherited
        
        self.model_type = model_type
            
    def run_iteration(self):
        '''
        Run a single deterministic OpenDrift run as part of a stochasitic ensemble
        This includes setting up the run directory and using the template config.py
        file to create a new one with an updated release_start_time
        '''
        # ---------------------------------------------------
        # prepare the iteration directory and config.py file
        # ---------------------------------------------------
        #
        os.makedirs(self.iteration_dir, exist_ok=True)
        os.chdir(self.iteration_dir)
        #
        config_template = os.path.join(self.run_dir,'config.py')
        config_destination = os.path.join(self.iteration_dir,'config.py')
        line_releasetime = 'release_start_time=\''+self.iteration_dir[-12:-1]+'\''
        # Open the template file for reading and create the new file for writing
        with open(config_template, 'r') as infile, open(config_destination, 'w') as outfile:
          for line in infile:
            if "release_start_time" in line:
              outfile.write(f"{line_releasetime}\n")
            else:
              outfile.write(line)
        
        # --------------
        # Run the model
        # --------------
        #
        if self.model_type == 'oil':
            run_oil(self.iteration_dir)
        elif self.model_type == 'leeway':
            run_leeway(self.iteration_dir)
        elif self.model_type == 'oceandrift':
            run_oceandrift(self.iteration_dir)
        # elif other model types...
        else:
            print('model_type not recognised: ' + self.model_type)
    
    def run_all(self):
        '''
        Loop through the specified run_id's and run the model for each'
        '''
        while self.run_id <= self.run_id_end:
            print('working on run'+ str(self.run_id).zfill(3))
            self.run_iteration()
            self.run_id += 1
            self.update_iteration_dir()

class grid_stochastic(base_stochastic):
    '''
    class used to go from particle OpenDrift output to gridded concentrations as part of a stochastic analysis
    '''
    def __init__(self, run_dir, date_start, run_id, increment_days, run_id_end,
                 fname_gridded='gridded.nc',
                 grid_type='density',
                 extents=None,
                 dx_m=None,
                 max_only=True):
        super().__init__(run_dir, date_start, run_id, increment_days, run_id_end) # this is where the base_stochastic class is inherited
        
        self.fname_gridded = fname_gridded
        self.grid_type        = grid_type
        self.extents          = extents
        self.dx_m             = dx_m
        self.max_only         = max_only
            
    def grid_iteration(self):
        '''
        Grid the output of a single deterministic OpenDrift run as part of a stochasitic ensemble
        '''
        fname = os.path.join(self.iteration_dir,'trajectories.nc')
        fname_gridded = os.path.join(self.iteration_dir,self.fname_gridded)
        grid_particles(fname,fname_gridded,grid_type=self.grid_type,extents=self.extents,dx_m=self.dx_m,max_only=self.max_only)
    
    def grid_all(self):
        '''
        Loop through the specified run_id's and run the model for each'
        '''
        while self.run_id <= self.run_id_end:
            self.grid_iteration()
            self.run_id += 1
            self.update_iteration_dir()

class gridded_stats(base_stochastic):
    '''
    class used to compute the probability of occurrence over a defined threshold and minimum time to arrival
    (assuming that the trajectories files have already been gridded)
    '''
    def __init__(self, run_dir, date_start, run_id, increment_days, run_id_end, 
                 out_dir = '/summary_stats/', # this gets appended onto run_dir
                 fname_gridded='gridded.nc', # the gridded filename common to all run directories
                 threshold=0 # only data over this value are used in computing statistics
                 ):
        super().__init__(run_dir, date_start, run_id, increment_days, run_id_end) # this is where the base_stochastic class is inherited
      
        self.out_dir = os.path.join(run_dir,out_dir)
        os.makedirs(self.out_dir, exist_ok=True)
        self.threshold              = threshold # only grid cells which exceed this value are counted
        self.fname_gridded       = fname_gridded # the generic gridded file name
        
        # get the grid from the first file
        with xr.open_dataset(self.iteration_dir+self.fname_gridded) as ds:
            lon_bin = ds.lon_bin.data
            lat_bin = ds.lat_bin.data
            
        # initialise an xarray dataset used for counting exceedance of threshold
        data_zeros=np.zeros((len(lat_bin), len(lon_bin)))
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
    
    def update_stats(self):
        '''
        read the output from one of the deterministic runs and update the probability and minimum time variables
        '''
        # get the data for this iteration
        with xr.open_dataset(self.iteration_dir+self.fname_gridded) as ds:
            iteration_max = ds.maximum.data # 'maximum' is the standard variable name accross all gridded files, representing the maximum over the deterministic simulation
            iteration_mintime =ds.minimum_time.data
        # update maximum if data from this iteration is greater
        previous_max=self.stats.maximum.data
        self.stats.maximum.data=np.maximum(iteration_max,previous_max)
        # save the run_id where this iteration is the overall maximum
        run_id_max=self.stats.run_id_of_max.data
        run_id_max[iteration_max>previous_max]=self.run_id
        self.stats.run_id_of_max.data=run_id_max
        # check where iteration max exceeds threshold and assign a value of 1
        iteration_exceed = np.zeros_like(iteration_max)
        iteration_exceed[iteration_max>self.threshold]=1
        # add up all the grid cells where the threshold is exceeded
        if len(self.no_grid_cells_over_threshold)==0:
            self.no_grid_cells_over_threshold=np.atleast_1d(np.sum(iteration_exceed))
            self.run_ids=np.atleast_1d(self.run_id)
        else:
            self.no_grid_cells_over_threshold=np.concatenate((self.no_grid_cells_over_threshold,np.atleast_1d(np.sum(iteration_exceed))),axis=0)
            self.run_ids=np.concatenate((self.run_ids,np.atleast_1d(self.run_id)),axis=0)
        self.stats.attrs["no_grid_cells_over_threshold"] = self.no_grid_cells_over_threshold # so that it is accessible in the output netcdf file we save later
        self.stats.attrs["run_ids"] =  self.run_ids
        # update count where threshold is exceeded
        self.stats.count_exceed.data += iteration_exceed
        # update probability of exceedance
        prob=self.stats.count_exceed.data/self.num_it
        # prob[prob==0]=np.nan
        self.stats.probability.data = prob
        
        # update the minimum time
        iteration_minimum_time = xr.where(iteration_exceed==1,iteration_mintime,10e6) 
        # update minimum time if data from this iteration is smaller
        previous_minimum_time=self.stats.minimum_time.data
        self.stats.minimum_time.data=np.minimum(iteration_minimum_time,previous_minimum_time)
    
    def update_stats_all(self):
        '''
        Loop through the specified run_id's and compute the stats over all runs
        '''
        while self.run_id <= self.run_id_end:
            self.update_stats()
            self.run_id += 1
            self.update_iteration_dir()
            self.num_it += 1
        self.__to_netcdf__()
    
    def __to_netcdf__(self):
        fname_out = self.fname_gridded.split('.nc')[0]+'_threshold'+str(self.threshold)+'.nc'
        self.stats.to_netcdf(os.path.join(self.out_dir,fname_out))

class gridded_stats_polygon(base_stochastic):
    '''
    class used to compute the probability of occurrence over a defined threshold and minimum time to arrival
    within a particular polygon
    (assuming that the trajectories files have already been gridded 
     - only gridded data inside the polygon are considered)
    '''
    def __init__(self, run_dir, date_start, run_id, increment_days, run_id_end, 
                 out_dir = '/summary_stats/', # this gets appended onto run_dir
                 fname_gridded='gridded.nc', # the gridded filename common to all run directories
                 threshold=0, # only data over this value are used in computing statistics
                 polygon_file='polygon_file.txt', # with 2 columns - lon lat
                 fname_out='polygon_stats.txt' # the output file
                 ):
        super().__init__(run_dir, date_start, run_id, increment_days, run_id_end) # this is where the base_stochastic class is inherited
      
        self.out_dir = os.path.join(run_dir,out_dir)
        os.makedirs(self.out_dir, exist_ok=True)
        self.threshold     = threshold # only grid cells which exceed this value are counted
        self.fname_gridded = fname_gridded # the generic gridded file name
        self.polygon_file  = polygon_file
        self.fname_out     = fname_out
        
        # Create mask for points inside the polygon
        with xr.open_dataset(self.iteration_dir+self.fname_gridded) as ds:
            lon = ds.lon_bin.data
            lat = ds.lat_bin.data
            lon_grid, lat_grid = np.meshgrid(lon, lat)
            lon = lon_grid.ravel()
            lat = lat_grid.ravel()

        polygon = np.loadtxt(polygon_file)
        # Define polygon path for masking
        polygon_path = Path(polygon)
        points = np.vstack((lon, lat)).T
        self.inside_polygon_mask = polygon_path.contains_points(points)
        
        # initialise variables used for counting exceedance of threshold
        self.count_exceed=0
        self.probability=0
        self.maximum=0
        self.minimum_time=10e6
    
    def update_stats(self):
        '''
        read the output from one of the deterministic runs and update the probability and minimum time variables
        '''
        # get the data for this iteration, inside the polygon
        with xr.open_dataset(self.iteration_dir+self.fname_gridded) as ds:     
            iteration_max = ds.maximum.data.ravel() # 'maximum' is the standard variable name accross all gridded files, representing the maximum over the deterministic simulation
            iteration_mintime =ds.minimum_time.data.ravel()
            # extract the max/min_time from points inside the polygon
            iteration_max=np.nanmax(iteration_max[self.inside_polygon_mask])
            iteration_mintime=np.nanmin(iteration_mintime[self.inside_polygon_mask])
            
        # update maximum if data from this iteration is greater
        self.maximum=np.maximum(iteration_max,self.maximum)
        # check if iteration max exceeds threshold and update the probability
        if iteration_max>self.threshold:
            self.count_exceed += 1
            self.minimum_time=np.minimum(iteration_mintime,self.minimum_time)
        # update probability of exceedance
        prob=self.count_exceed/self.num_it
        self.probability = prob
    
    def update_stats_all(self):
        '''
        Loop through the specified run_id's and compute the stats over all runs
        '''
        while self.run_id <= self.run_id_end:
            self.update_stats()
            self.run_id += 1
            self.update_iteration_dir()
            self.num_it += 1
        
        # write the stats to the output txt file
        fname_out = os.path.join(self.out_dir,self.fname_out)
        with open(fname_out, 'w') as f:
            f.write(f"probability={self.probability}\n")
            f.write(f"maximum={self.maximum}\n")
            f.write(f"minimum_time={self.minimum_time}\n")    
    
class stochasitic_massbal(base_stochastic):
    '''
    class used to compute the stochastic mass balance
    '''
    def __init__(self, run_dir, date_start, run_id, increment_days, run_id_end, 
                 out_dir = '/summary_stats/', # this gets appended onto run_dir
                 fname_massbal='oil_massbal.nc'):
        super().__init__(run_dir, date_start, run_id, increment_days, run_id_end) # this is where the base_stochastic class is inherited
        
        self.out_dir = os.path.join(run_dir,out_dir)
        os.makedirs(self.out_dir, exist_ok=True)
        self.fname_massbal = fname_massbal
        # initialise with the mass balance of the first run_id
        self.massbal_iteration()
        ds = self.get_ds_massbal()
        self.budget = ds
    
    def massbal_iteration(self):
        '''
        do the mass balance of a single deterministic OpenOil run as part of a stochasitic ensemble
        '''
        fname = os.path.join(self.iteration_dir,'trajectories.nc')
        fname_massbal = os.path.join(self.iteration_dir,self.fname_massbal)
        oil_massbal(fname,fname_massbal)
    
    def get_ds_massbal(self):
        with xr.open_dataset(os.path.join(self.iteration_dir,self.fname_massbal)) as ds:
            # replace the times with time since start of run for comparison with other runs
            days_since_start = (ds.time.data-ds.time.data[0]).astype('timedelta64[s]').astype(np.int32)/3600/24 # seriously convoluted but it works. Timedelta64 is not ideals
            ds = ds.assign(time = days_since_start)
            return ds
    
    def update_massbal_all(self):
        '''
        Loop through the specified run_id's and compute the mass balance over all runs
        '''
        while self.run_id <= self.run_id_end:
            if self.num_it>1: # but only if we're not on the first iteration (this has already been added in __init__)
                self.massbal_iteration()
                ds = self.get_ds_massbal()
                self.budget=xr.concat([self.budget, ds], dim="iteration")
            self.run_id += 1
            self.update_iteration_dir()
            self.num_it += 1
        # fill missing values with the last valid data
        # (to account for simulations which finish early due to all particles being deactivated)
        self.budget=self.budget.ffill('time')
        self.__to_netcdf__()
        
    def __to_netcdf__(self):
        self.budget.to_netcdf(os.path.join(self.out_dir,self.fname_massbal))
        
#if __name__ == "__main__":
    
    
