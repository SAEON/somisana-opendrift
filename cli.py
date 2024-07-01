'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow python functions to be run from the cli docker image for this repo
So this is the entry point for the docker image (see the Dockerfile).
But it's also handy if you want to execute a python function from inside a bash script
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import sys
import argparse
from opendrift_tools.run import oil as run_oil
from opendrift_tools.postprocess import grid_particles
from opendrift_tools.plotting import plot_particles
from opendrift_tools.plotting import plot_gridded

def parse_list(value):
    return [x.strip() for x in value.split(',')]

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-opendrift repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...

    # ----------------
    # run_oil
    # ----------------
    parser_run_oil = subparsers.add_parser('run_oil', 
            help='Run an OpenOil simulation')
    parser_run_oil.add_argument('--config_dir', required=True, type=str, help='Directory where the config.py file is located')
    def run_oil_handler(args):
        run_oil(args.config_dir)
    parser_run_oil.set_defaults(func=run_oil_handler)
    
    # -------------------------
    # grid the particle output
    # -------------------------
    parser_grid_particles = subparsers.add_parser('grid_particles', 
            help='convert the particle output of an OpenDrift simulation to a eulerian grid')
    parser_grid_particles.add_argument('--config_dir', required=True, type=str, help='Directory where the config.py file is located')
    def grid_particles_handler(args):
        sys.path.append(args.config_dir)
        import config
        grid_particles(config.fname,
                       config.fname_gridded,
                       extents=config.grid_extents,
                       dx_m=config.dx_m,
                       max_only=config.max_only)
    parser_grid_particles.set_defaults(func=grid_particles_handler)
    
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
        plot_particles(config.fname,
                        figsize=config.figsize,
                        extents=config.plot_extents,
                        lscale=config.lscale,
                        lon_release=config.lon_release,
                        lat_release=config.lat_release,
                        size_release=config.size_release,
                        size_scat=config.size_scat,
                        gif_out=config.gif_out_particles,
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
        plot_gridded(config.fname_gridded,
                        figsize=config.figsize,
                        extents=config.plot_extents,
                        lscale=config.lscale,
                        lon_release=config.lon_release,
                        lat_release=config.lat_release,
                        size_release=config.size_release,
                        gif_out=config.gif_out_gridded,
                        write_gif=config.write_gif,
                        skip_time=config.skip_time,
                        tstep_end=config.tstep_end)
    parser_plot_gridded.set_defaults(func=plot_gridded_handler)
    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    
