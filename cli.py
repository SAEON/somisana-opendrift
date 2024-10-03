'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow python functions to be run from the cli docker image for this repo
So this is the entry point for the docker image (see the Dockerfile).
But it's also handy if you want to execute a python function from inside a bash script
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import sys, os
import argparse
from datetime import datetime
from opendrift_tools.run import oil as run_oil
from opendrift_tools.run import leeway as run_leeway
from opendrift_tools.run import oceandrift as run_oceandrift
from opendrift_tools.postprocess import grid_particles, oil_massbal, combine_gridded, combine_trajectories
from opendrift_tools.plotting import plot_particles, plot_gridded, plot_gridded_stats, plot_budget
from opendrift_tools.stochastic import run_stochastic, grid_stochastic, gridded_stats, gridded_stats_polygon, stochasitic_massbal

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y%m%d_%H')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYYMMDD_HH'.")

def parse_list(value):
    if value is None or value == 'None':
        return None
    else:
        return [float(x.strip()) for x in value.split(',')]

def parse_list_str(value):
    if value is None or value == 'None':
        return None
    else:
        return [x.strip() for x in value.split(',')]

def parse_float(value):
    if value is None or value == 'None':
        return None
    else:
        return float(value)

def parse_bool(s: str) -> bool:
    try:
        return {'true':True, 'false':False}[s.lower()]
    except KeyError:
        raise argparse.ArgumentTypeError(f'expect true/false, got: {s}')

def parse_str(value):
    if value is None:
        return None
    try:
        return str(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid string: {value}")

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-opendrift repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...

    # ----------------
    # run_model
    # ----------------
    parser_run_model = subparsers.add_parser('run_model', 
            help='Run an OpenDrift simulation')
    parser_run_model.add_argument('--config_dir', required=True, type=str, help='Directory where the config.py file is located')
    parser_run_model.add_argument('--model_type', required=True, type=str, default='oil', help='type of model to run- options are \'oceandrift\', \'oil\' or \'leeway\'')
    def run_model_handler(args):
        if args.model_type == 'oil':
            run_oil(args.config_dir)
        elif args.model_type == 'leeway':
            run_leeway(args.config_dir)
        elif args.model_type == 'oceandrift':
            run_oceandrift(args.config_dir)
        else:
            raise argparse.ArgumentTypeError(f"Invalid model_type: {args.model_type}")
    parser_run_model.set_defaults(func=run_model_handler)
    
    # -------------------------
    # grid the particle output
    # -------------------------
    parser_grid_particles = subparsers.add_parser('grid_particles', 
            help='convert the particle output of an OpenDrift simulation to a eulerian grid')
    parser_grid_particles.add_argument('--config_dir', required=True, type=str, help='Directory where the OpenDrift output is located')
    parser_grid_particles.add_argument('--fname', required=False, type=str, default='trajectories.nc', help='the OpenDrift output filename')
    parser_grid_particles.add_argument('--fname_gridded', required=False, type=str, default='gridded.nc', help='the gridded filename')
    parser_grid_particles.add_argument('--grid_type', required=False, type=str, default='density', help='what kind of gridding to do. Options are \'density\', \'surface_oil\', \'stranded_oil\'')
    parser_grid_particles.add_argument('--extents', required=False,type=parse_list, default=None, help='the spatial extent of the grid in format lon0,lon1,lat0,lat1. If None, then this is automatically determined from the geographic extent of the particles')
    parser_grid_particles.add_argument('--dx_m', required=False, type=parse_float, default=None, help='grid size in meters. If None, then a 50 x 50 regular grid is generated')
    parser_grid_particles.add_argument('--max_only', required=False, type=parse_bool, default=False,
            help='option to only write the maximum over the entire file to save disk space (set to true or false)')
    def grid_particles_handler(args):
        fname = os.path.join(args.config_dir,args.fname)
        fname_gridded = os.path.join(args.config_dir,args.fname_gridded)
        grid_particles(fname,
                       fname_gridded,
                       grid_type=args.grid_type,
                       extents=args.extents,
                       dx_m=args.dx_m,
                       max_only=args.max_only)
    parser_grid_particles.set_defaults(func=grid_particles_handler)
    
    # ----------------------------------------------------
    # combine multiple gridded outputs into a single file
    # ----------------------------------------------------
    parser_combine_gridded = subparsers.add_parser('combine_gridded', 
            help='combine gridded output from multiple simulations into a single file')
    parser_combine_gridded.add_argument('--dir_with_dirs', required=True, type=str, help='path to a directory containing multiple directories, where each directory contains opendrift output')
    parser_combine_gridded.add_argument('--dirs', required=False, type=parse_list_str, default=None, help='a comma separated list of dir names inside dir_with_dirs which contain opendrift ouput. If None, then dirs is found dynamically from all dirs inside dir_with_dirs (except for the dirname \'combined\')')
    parser_combine_gridded.add_argument('--fname_gridded', required=False, type=str, default='gridded.nc', help='the gridded filename to be found in every dir in dirs')
    parser_combine_gridded.add_argument('--var_name', required=False, type=parse_str, default=None, help='the time-dependent variable in fname_gridded. If None, then only the time-independent \'maximum\' and \'minimum_time\' variables are used')
    def combine_gridded_handler(args):
        combine_gridded(args.dir_with_dirs,
                       dirs=args.dirs,
                       fname_gridded=args.fname_gridded,
                       var_name=args.var_name)
    parser_combine_gridded.set_defaults(func=combine_gridded_handler)
    
    # ---------------------------------------------------------
    # combine multiple trajectories outputs into a single file
    # ---------------------------------------------------------
    parser_combine_trajectories = subparsers.add_parser('combine_trajectories', 
            help='combine raw opendrift output from multiple simulations into a single file')
    parser_combine_trajectories.add_argument('--dir_with_dirs', required=True, type=str, help='path to a directory containing multiple directories, where each directory contains opendrift output')
    parser_combine_trajectories.add_argument('--dirs', required=False, type=parse_list_str, default=None, help='a comma separated list of dir names inside dir_with_dirs which contain opendrift ouput. If None, then dirs is found dynamically from all dirs inside dir_with_dirs (except for the dirname \'combined\')')
    parser_combine_trajectories.add_argument('--fname_traj', required=False, type=str, default='trajectories.nc', help='the output opendrift filename to be found in every dir in dirs')
    def combine_trajectories_handler(args):
        combine_trajectories(args.dir_with_dirs,
                       dirs=args.dirs,
                       fname_traj=args.fname_traj)
    parser_combine_trajectories.set_defaults(func=combine_trajectories_handler)
    
    # -------------------------------------------
    # get the oil mass balance of an OpenOil run
    # -------------------------------------------
    parser_oil_massbal = subparsers.add_parser('oil_massbal', 
            help='compute the mass balance for an OpenOil simulation')
    parser_oil_massbal.add_argument('--config_dir', required=True, type=str, help='Directory where the OpenOil output is located')
    parser_oil_massbal.add_argument('--fname', required=False, type=str, default='trajectories.nc', help='the OpenOil output filename')
    parser_oil_massbal.add_argument('--fname_out', required=False, type=str, default='oil_mass_balance.nc', help='the mass balance filename')
    def oil_massbal_handler(args):
        fname = os.path.join(args.config_dir,args.fname)
        fname_out = os.path.join(args.config_dir,args.fname_out)
        oil_massbal(fname, fname_out)
    parser_oil_massbal.set_defaults(func=oil_massbal_handler)
    
    # -------------------------------------------
    # plot the oil mass balance of an OpenOil run
    # -------------------------------------------
    parser_plot_budget = subparsers.add_parser('plot_budget', 
            help='plot the mass balance for an OpenOil simulation, as output from the oil_massbal function')
    parser_plot_budget.add_argument('--config_dir', required=True, type=str, help='Directory where the OpenOil output is located')
    parser_plot_budget.add_argument('--fname', required=False, type=str, default='oil_mass_balance.nc', help='the mass balance filename output from the oil_massbal function')
    parser_plot_budget.add_argument('--fname_out', required=False, type=str, default='oil_mass_balance.jpg', help='the output plot filename')
    def plot_budget_handler(args):
        fname = os.path.join(args.config_dir,args.fname)
        fname_out = os.path.join(args.config_dir,args.fname_out)
        plot_budget(fname, fname_out)
    parser_plot_budget.set_defaults(func=plot_budget_handler)
    
    # ---------------------------------------------------
    # do a plot/animation of the particle/gridded output
    # ---------------------------------------------------
    parser_animate = subparsers.add_parser('animate', 
            help='do a plot/animation of the particle/gridded output')
    parser_animate.add_argument('--type', required=True, type=str, help='either particles, gridded or gridded_stats')
    parser_animate.add_argument('--config_dir', required=True, type=str, help='Directory where the output files are located')
    parser_animate.add_argument('--fname_gridded', required=False, type=str, default='gridded.nc', help='the gridded filename')
    parser_animate.add_argument('--fname_particles', required=False, type=str, default='trajectories.nc', help='the raw OpenDrift particle output filename')
    parser_animate.add_argument('--var_str', required=False, default='particle_density', type=str, help='the variable name to plot')
    parser_animate.add_argument('--extents', required=False,type=parse_list, default=None, help='the spatial extent of the plot in format lon0,lon1,lat0,lat1. If None, then this is automatically determined from the geographic extent of the particles')
    parser_animate.add_argument('--lscale', required=False, type=str, default='h', help='land resolution, with options c,l,i,h,f and auto')
    parser_animate.add_argument('--ticks', required=False, type=parse_list,
                         default=[0,1,2,3,5,10,15],
                         help='contour ticks to use in plotting the variable')
    parser_animate.add_argument('--cbar_label', required=False, default='particle density (%)', type=str, help='the label used for the colorbar')
    parser_animate.add_argument('--cmap', required=False, default='Spectral_r', type=str, help='the colormap to use')
    parser_animate.add_argument('--gif_out', required=False, type=parse_str, default=None, help='the output gif filename')
    parser_animate.add_argument('--jpg_out', required=False, type=parse_str, default=None, help='the output jpg filename')
    def animate_handler(args):
        # the input options passed by the cli is not exhaustive to the plotting functions
        # this is just intended to provide a quick animation as part of the operational workflow
        # the input cli arguments could easily be extended if deemed necessary
        fname_gridded = os.path.join(args.config_dir,args.fname_gridded)
        fname_particles = os.path.join(args.config_dir,args.fname_particles)
        if args.gif_out is None: # we need this if/else to ensure the os.path.join() gets a string, not a None
            gif_out = None
        else:
            gif_out = os.path.join(args.config_dir,args.gif_out)
        if args.jpg_out is None: # we need this if/else to ensure the os.path.join() gets a string, not a None
            jpg_out = None
        else:
            jpg_out = os.path.join(args.config_dir,args.jpg_out)
        
        if args.type == 'particles':
            # animate the particles
            plot_particles(fname_particles,
                            extents=args.extents,
                            lscale=args.lscale,
                            gif_out=gif_out,
                            write_gif=True)
        elif args.type == 'gridded':
            # animate the gridded output
            plot_gridded(fname_gridded,
                            fname_particles=fname_particles,
                            var_str=args.var_str,
                            extents=args.extents,
                            lscale=args.lscale,
                            ticks=args.ticks,
                            cbar_label=args.cbar_label,
                            cmap=args.cmap,
                            gif_out=gif_out,
                            write_gif=True)
        elif args.type == 'gridded_stats':
            # plot the static output (i.e. no time dimension) of the gridded output
            plot_gridded_stats(fname_gridded,
                            fname_particles=fname_particles,
                            var_str=args.var_str,
                            extents=args.extents,
                            lscale=args.lscale,
                            ticks=args.ticks,
                            cbar_label=args.cbar_label,
                            cmap=args.cmap,
                            jpg_out=jpg_out,
                            write_jpg=True)
        else:
            raise argparse.ArgumentTypeError(f"Invalid animation type: {args.type}")
    parser_animate.set_defaults(func=animate_handler)
    
    # ----------------
    # run_stochastic
    # ----------------
    parser_run_stochastic = subparsers.add_parser('run_stochastic', 
            help='Run stochastic OpenDrift simulations')
    parser_run_stochastic.add_argument('--run_dir', required=True, type=str, help='based dir where stochastic iterations are initialised')
    parser_run_stochastic.add_argument('--date_start', required=True, type=parse_datetime, help='start time of run001 (first stochastic simulation) in format "YYYYMMDD_HH"')
    parser_run_stochastic.add_argument('--run_id', required=True, type=int, help='run id to start on (you don\'t have to start at run001)')
    parser_run_stochastic.add_argument('--increment_days', required=True, type=float, help='number of days increment between stochastic runs')
    parser_run_stochastic.add_argument('--run_id_end', required=True, type=int, help='run id to end on')
    parser_run_stochastic.add_argument('--model_type', required=False, type=str, default='oil', help='type of model to run- options are \'oceandrift\', \'oil\' or \'leeway\'')
    def run_stochastic_handler(args):
        stoch = run_stochastic(args.run_dir, args.date_start, args.run_id, args.increment_days, args.run_id_end, model_type=args.model_type)
        stoch.run_all()
    parser_run_stochastic.set_defaults(func=run_stochastic_handler)
    
    # ----------------
    # grid_stochastic
    # ----------------
    parser_grid_stochastic = subparsers.add_parser('grid_stochastic', 
            help='Grid the output of stochastic OpenDrift simulations')
    parser_grid_stochastic.add_argument('--run_dir', required=True, type=str, help='based dir where stochastic iterations are initialised')
    parser_grid_stochastic.add_argument('--date_start', required=True, type=parse_datetime, help='start time of run001 (first stochastic simulation) in format "YYYYMMDD_HH"')
    parser_grid_stochastic.add_argument('--run_id', required=True, type=int, help='run id to start on (you don\'t have to start at run001)')
    parser_grid_stochastic.add_argument('--increment_days', required=True, type=float, help='number of days increment between stochastic runs')
    parser_grid_stochastic.add_argument('--run_id_end', required=True, type=int, help='run id to end on')
    parser_grid_stochastic.add_argument('--fname_gridded', required=False, type=str, default='gridded.nc', help='the gridded filename')
    parser_grid_stochastic.add_argument('--grid_type', required=False, type=str, default='density', help='what kind of gridding to do. Options are \'density\', \'surface_oil\', \'stranded_oil\'')
    parser_grid_stochastic.add_argument('--extents', required=False,type=parse_list, default=None, help='the spatial extent of the grid in format lon0,lon1,lat0,lat1. If None, then this is automatically determined from the geographic extent of the particles')
    parser_grid_stochastic.add_argument('--dx_m', required=False, type=float, default=None, help='grid size in meters. If None, then a 50 x 50 regular grid is generated')
    parser_grid_stochastic.add_argument('--max_only', required=False, type=parse_bool, default=False,
            help='option to only write the maximum over the entire file to save disk space (set to true or false)')
    def grid_stochastic_handler(args):
        stoch = grid_stochastic(args.run_dir, args.date_start, args.run_id, args.increment_days, args.run_id_end, 
                                fname_gridded=args.fname_gridded,
                                grid_type=args.grid_type,
                                extents=args.extents,
                                dx_m=args.dx_m,
                                max_only=args.max_only)
        stoch.grid_all()
    parser_grid_stochastic.set_defaults(func=grid_stochastic_handler)
    
    # --------------
    # gridded_stats
    # --------------
    parser_gridded_stats = subparsers.add_parser('gridded_stats', 
            help='Compute statistics on gridded output from stochastic OpenDrift simulations')
    parser_gridded_stats.add_argument('--run_dir', required=True, type=str, help='based dir where stochastic iterations are initialised')
    parser_gridded_stats.add_argument('--date_start', required=True, type=parse_datetime, help='start time of run001 (first stochastic simulation) in format "YYYYMMDD_HH"')
    parser_gridded_stats.add_argument('--run_id', required=True, type=int, help='run id to start on (you don\'t have to start at run001)')
    parser_gridded_stats.add_argument('--increment_days', required=True, type=float, help='number of days increment between stochastic runs')
    parser_gridded_stats.add_argument('--run_id_end', required=True, type=int, help='run id to end on')
    parser_gridded_stats.add_argument('--out_dir', required=False, type=str,default='summary_stats', help='output directory name (this gets appended onto run_dir)')
    parser_gridded_stats.add_argument('--fname_gridded', required=True, type=str, help='the gridded filename common to all run directories')
    parser_gridded_stats.add_argument('--threshold', required=True, type=float, help='only data over this value are used in computing statistics')
    def gridded_stats_handler(args):
        stoch = gridded_stats(args.run_dir, args.date_start, args.run_id, args.increment_days, args.run_id_end, 
                                out_dir = args.out_dir,
                                fname_gridded=args.fname_gridded,
                                threshold=args.threshold)
        stoch.update_stats_all()
    parser_gridded_stats.set_defaults(func=gridded_stats_handler)
    
    # ----------------------
    # gridded_stats_polygon
    # ----------------------
    parser_gridded_stats_polygon = subparsers.add_parser('gridded_stats_polygon', 
            help='Compute summary statistics within a polygon from gridded output from stochastic OpenDrift simulations')
    parser_gridded_stats_polygon.add_argument('--run_dir', required=True, type=str, help='based dir where stochastic iterations are initialised')
    parser_gridded_stats_polygon.add_argument('--date_start', required=True, type=parse_datetime, help='start time of run001 (first stochastic simulation) in format "YYYYMMDD_HH"')
    parser_gridded_stats_polygon.add_argument('--run_id', required=True, type=int, help='run id to start on (you don\'t have to start at run001)')
    parser_gridded_stats_polygon.add_argument('--increment_days', required=True, type=float, help='number of days increment between stochastic runs')
    parser_gridded_stats_polygon.add_argument('--run_id_end', required=True, type=int, help='run id to end on')
    parser_gridded_stats_polygon.add_argument('--out_dir', required=False, type=str,default='summary_stats', help='output directory name (this gets appended onto run_dir)')
    parser_gridded_stats_polygon.add_argument('--fname_gridded', required=True, type=str, help='the gridded filename common to all run directories')
    parser_gridded_stats_polygon.add_argument('--threshold', required=True, type=float, help='only data over this value are used in computing statistics')
    parser_gridded_stats_polygon.add_argument('--polygon_file', required=True, type=str, help='full path to the input polygon text file. First column is lon, second column is lat. Space delimited')
    parser_gridded_stats_polygon.add_argument('--fname_out', required=True, type=str, help='the output filename (text file) to be written to out_dir')
    def gridded_stats_polygon_handler(args):
        stoch = gridded_stats_polygon(args.run_dir, args.date_start, args.run_id, args.increment_days, args.run_id_end, 
                                out_dir = args.out_dir,
                                fname_gridded=args.fname_gridded,
                                threshold=args.threshold,
                                polygon_file=args.polygon_file,
                                fname_out=args.fname_out)
        stoch.update_stats_all()
    parser_gridded_stats_polygon.set_defaults(func=gridded_stats_polygon_handler)
    
    # -----------------------
    # stohastic mass balance
    # -----------------------
    parser_stoch_massbal = subparsers.add_parser('stochastic_massbal', 
            help='Compute the stochasitc mass balance from stochastic OpenOil simulations')
    parser_stoch_massbal.add_argument('--run_dir', required=True, type=str, help='based dir where stochastic iterations are initialised')
    parser_stoch_massbal.add_argument('--date_start', required=True, type=parse_datetime, help='start time of run001 (first stochastic simulation) in format "YYYYMMDD_HH"')
    parser_stoch_massbal.add_argument('--run_id', required=True, type=int, help='run id to start on (you don\'t have to start at run001)')
    parser_stoch_massbal.add_argument('--increment_days', required=True, type=float, help='number of days increment between stochastic runs')
    parser_stoch_massbal.add_argument('--run_id_end', required=True, type=int, help='run id to end on')
    parser_stoch_massbal.add_argument('--out_dir', required=False, type=str,default='summary_stats', help='output directory name (this gets appended onto run_dir)')
    parser_stoch_massbal.add_argument('--fname_massbal', required=False, type=str, default='oil_massbal.nc', help='the mass balance filename')
    def stoch_massbal_handler(args):
        stoch = stochasitic_massbal(args.run_dir, args.date_start, args.run_id, args.increment_days, args.run_id_end, 
                                out_dir = args.out_dir,
                                fname_massbal=args.fname_massbal
                                )
        stoch.update_massbal_all()
    parser_stoch_massbal.set_defaults(func=stoch_massbal_handler)
    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    
