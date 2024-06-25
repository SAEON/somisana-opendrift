'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow python functions to be run from the cli docker image for this repo
So this is the entry point for the docker image (see the Dockerfile).
But it's also handy if you want to execute a python function from inside a bash script
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import os
import argparse
from opendrift_tools.run import oil as run_oil

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
    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    
