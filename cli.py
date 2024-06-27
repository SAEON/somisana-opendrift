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
from opendrift_tools.plotting import plot_particles

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
    parser_run_oil.add_argument('--config_dir', required=True, type=str, help='Directory where the config_oil.py file is located')
    def run_oil_handler(args):
        run_oil(args.config_dir)
    parser_run_oil.set_defaults(func=run_oil_handler)
    
    # -----------------------
    # do a plot or animation
    # -----------------------
    parser_plot_particles = subparsers.add_parser('plot_particles', 
            help='do a plot or an animation of the particle output of an OpenDrift simulation')
    parser_plot_particles.add_argument('--config_dir', required=True, type=str, help='Directory where the config_oil.py file is located')
    def plot_particles_handler(args):
        # for now I would prefer the functionality of parsing all the info in the config file to the plot_particles() function
        # rather than just reading all the variables from the config file inside the function as it is done in run.py
        # maybe we should change this at some point and use the config file for everything to minimise function inputs?
        sys.path.append(args.config_dir)
        import config
        plot_particles(config.fname,
                        var_str=config.var_str,
                        tstep=config.tstep,
                        figsize=config.figsize,
                        extents=config.extents,
                        lscale=config.lscale,
                        lon_release=config.lon_release,
                        lat_release=config.lat_release,
                        size_release=config.size_release,
                        size_scat=config.size_scat,
                        ticks=config.ticks,
                        cmap=config.cmap,
                        plot_cbar=config.plot_cbar,
                        cbar_loc=config.cbar_loc,
                        cbar_label=config.cbar_label,
                        jpg_out=config.jpg_out,
                        write_jpg=config.write_jpg,
                        gif_out=config.gif_out,
                        write_gif=config.write_gif,
                        skip_time=config.skip_time,
                        tstep_end=config.tstep_end)
    parser_plot_particles.set_defaults(func=plot_particles_handler)
    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    
