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
from opendrift_tools.postprocess import grid_particles, oil_massbal
from opendrift_tools.plotting import plot_particles, plot_gridded
from opendrift_tools.stochastic import run_stochastic, grid_stochastic, gridded_stats, stochasitic_massbal

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y%m%d_%H')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYYMMDD_HH'.")

def parse_list(value):
    if value is not None:
        return [float(x.strip()) for x in value.split(',')]
    else:
        return None

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
            print('model_type not recognised: ' + args.model_type)
    parser_run_model.set_defaults(func=run_model_handler)
    
    # -------------------------
    # grid the particle output
    # -------------------------
    parser_grid_particles = subparsers.add_parser('grid_particles', 
            help='convert the particle output of an OpenDrift simulation to a eulerian grid')
    parser_grid_particles.add_argument('--config_dir', required=True, type=str, help='Directory where the OpenDrift output is located')
    parser_grid_particles.add_argument('--fname', required=False, type=str, default='trajectories.nc', help='the OpenDrift output filename')
    parser_grid_particles.add_argument('--fname_gridded', required=False, type=str, default='gridded_density.nc', help='the gridded filename')
    parser_grid_particles.add_argument('--grid_type', required=False, type=str, default='density', help='what kind of gridding to do. Options are \'density\', \'surface_oil\', \'stranded_oil\'')
    parser_grid_particles.add_argument('--extents', required=False,type=parse_list, default=None, help='the spatial extent of the grid in format lon0,lon1,lat0,lat1. If None, then this is automatically determined from the geographic extent of the particles')
    parser_grid_particles.add_argument('--dx_m', required=False, type=float, default=None, help='grid size in meters. If None, then a 50 x 50 regular grid is generated')
    parser_grid_particles.add_argument('--max_only', required=False, action='store_true', help='option to only write the maximum over the entire file to save disk space (set this flag to enable)')
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
    
    # -------------------------------------------
    # get the oil mass balance of an OpenOil run
    # -------------------------------------------
    parser_oil_massbal = subparsers.add_parser('oil_massbal', 
            help='compute the mass balance for an OpenOil simulation')
    parser_oil_massbal.add_argument('--config_dir', required=True, type=str, help='Directory where the OpenOil output is located')
    parser_oil_massbal.add_argument('--fname', required=False, type=str, default='trajectories.nc', help='the OpenOil output filename')
    parser_oil_massbal.add_argument('--fname_out', required=False, type=str, default='gridded.nc', help='the mass balance filename')
    def oil_massbal_handler(args):
        fname = os.path.join(args.config_dir,args.fname)
        fname_out = os.path.join(args.config_dir,args.fname_out)
        oil_massbal(fname, fname_out)
    parser_oil_massbal.set_defaults(func=oil_massbal_handler)
    
    # ----------------------------------------------
    # do a plot or animation of the particle output
    # ----------------------------------------------
    parser_plot_particles = subparsers.add_parser('plot_particles', 
            help='do a plot or an animation of the particle output of an OpenDrift simulation')
    parser_plot_particles.add_argument('--config_dir', required=True, type=str, help='Directory where the config.py file is located')
    def plot_particles_handler(args):
        # the input options passed by the cli is not exhaustive
        # just intended to provide a quick animation as part of the operational workflow
        sys.path.append(args.config_dir)
        import config
        fname = os.path.join(args.config_dir,config.fname)
        gif_out_particles = os.path.join(args.config_dir,config.gif_out_particles)
        plot_particles(fname,
                        figsize=config.figsize,
                        extents=config.plot_extents,
                        lscale=config.lscale,
                        lon_release=config.lon_release,
                        lat_release=config.lat_release,
                        size_release=config.size_release,
                        size_scat=config.size_scat,
                        gif_out=gif_out_particles,
                        write_gif=config.write_gif,
                        skip_time=config.skip_time,
                        tstep_end=config.tstep_end)
    parser_plot_particles.set_defaults(func=plot_particles_handler)
    
    # ----------------------------------------------
    # do a plot or animation of the gridded output
    # ----------------------------------------------
    parser_plot_gridded = subparsers.add_parser('plot_gridded', 
            help='do a plot or an animation of the gridded output from the grid_particles function')
    parser_plot_gridded.add_argument('--config_dir', required=True, type=str, help='Directory where the config.py file is located')
    def plot_gridded_handler(args):
        # the input options passed by the cli is not exhaustive
        # just intended to provide a quick animation as part of the operational workflow
        sys.path.append(args.config_dir)
        import config
        fname_gridded = os.path.join(args.config_dir,config.fname_gridded)
        gif_out_gridded = os.path.join(args.config_dir,config.gif_out_gridded)
        plot_gridded(fname_gridded,
                        figsize=config.figsize,
                        extents=config.plot_extents,
                        lscale=config.lscale,
                        lon_release=config.lon_release,
                        lat_release=config.lat_release,
                        size_release=config.size_release,
                        gif_out=gif_out_gridded,
                        write_gif=config.write_gif,
                        skip_time=config.skip_time,
                        tstep_end=config.tstep_end)
    parser_plot_gridded.set_defaults(func=plot_gridded_handler)
    
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
    parser_grid_stochastic.add_argument('--max_only', required=False, action='store_true', help='option to only write the maximum over the entire file to save disk space (set this flag to enable)')
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
    
