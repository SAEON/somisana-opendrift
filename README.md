# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used to run [OpenDrift](https://opendrift.github.io/) simulations on top of our [operational ocean modelling system](https://github.com/SAEON/somisana-croco). The [OpenDrift](https://opendrift.github.io/) documentation is very good, and one should start there if new to OpenDrift. The primary purpose of this repo is to provide some tools to run OpenDrift simulations operationally on a server set up for this purposes e.g. in response to search and rescue or oil spill emergencies. It is also possible to clone this repo locally, in which case it provides a way of running your own local simulations with minimal setup.  

Directories in the repository: 
- `opendrift_tools`:   python functions for some preprocessing, running, postprocessing and plotting OpenDrift simulations
- `configs`:           configuration files used for running OpenDrift simulations - these files are used to define the variables used in various pre and postprocessing functions
- `.github/workflows`: github workflows for running OpenDrift operationally on a server set up for this purpose on DFFE's Marine Information Management System (MIMS) 

# Setting up the python environment for local use and/or development

To run OpenDrift simulations locally, you need to first have the `opendrift` mamba/conda environment set up, according to the [OpenDrift install instructions](https://opendrift.github.io/install.html).

Once you've set up your `opendrift` environment, you can clone this repo to your local machine (do this wherever you would like the code):

`git clone git@github.com:SAEON/somisana-opendrift.git`

Then navigate to the root directory of the repo, and install the code into your `opendrift` environment:

```sh
cd somisana-opendrift
conda activate opendrift
pip install --no-deps -e .
```

Now you'll be able to import all the functions inside the `opendrift_tools` directory without having to manually add this directory to your python path.

# Example of running a simulation locally

To run a local simulation from this repo, all that is needed is to set up a relevant `config_*.py` file which are in the `configs` directory. These default configuration files are set up to run an operational simulation on MIMS. So to run locally, you would need to copy the relevant default file to wherever you want to run the model, and then edit it accordingly - particularly the paths to the forcing files which will obviously not be relevant for your local run. For example, to run an oil spill simulation:

```sh
cp <path-to-somisana-opendrift>/configs/config_oil.py <your-config-dir>/config.py
# then edit the file as you like
```

To run the model, you can put this in your own python file and run it:

```sh
from opendrift_tools.run import oil as run_oil
run_oil(<your-config-dir>)
```

Or you can run the model from the command line interface (cli.py ) for this repo. The cli.py is really intended to provide the entry point to the docker image for this repo, so that we can get access to relevant functions without the need to set up an environment wherever we want to run the model.

```sh
python /home/gfearon/code/somisana-opendrift/cli.py run_model --model_type oil --config_dir <your-config-dir>
```

# gridding the output

It is often desirable to go from the Lagrangian particle output to a Eulerian field which shows the particle density, or maybe surface oil thickness in the case of an oil spill. So far, this repo contains functions for calculating the surface oil thickness, stranded oil concentration, or the particle density (which would be applicable to all types of runs, not just oil spills). To call the function from your own python script:

```sh
from opendrift_tools.postprocess import grid_particles
grid_particles(...) # check out the documentation in the function for the expected inputs
```

Or again, you can do this from the command line interface (cli.py), in which case the inputs to the function are read from the config.py file:

```sh
python /home/gfearon/code/somisana-opendrift/cli.py grid_particles --config_dir <your-config-dir> # there are a bunch of optional inputs to this function you can check out
```

# plotting the output

There is a basic function for doing a plot or an animation of your raw output. You can edit the relevant section in your config.py file, and then to do the plot/animation:

```sh
python /home/gfearon/code/somisana-opendrift/cli.py plot_particles --config_dir <your-config-dir>
```

Or to plot the gridded output:

```sh
python /home/gfearon/code/somisana-opendrift/cli.py plot_gridded --config_dir <your-config-dir>
```

But if you want to personalise your plot more, it's best to do your own thing, perhaps based on the `plot_particles()` and `plot_gridded()` functions in `opendrift_tools/plotting.py`... 

# Running the docker image

As mentioned above, the command line interface (cli.py) in the root directory of this repo allows one command line access to various functions in this repo. This is the entry point for the Docker image which can be built from this repo. Have a look at the `Dockerfile` to see how it is built. At the time of writing, we are building on top of the `ghcr.io/saeon/opendrift_v1.11.0` base image, which is created using Dockerfile.base in this repo (see the explanatory text in there).

To build the docker image locally you could do:

```sh
docker build -t somisana_opendrift .
```

This creates a local docker image which has a default behaviour of running `cli.py`. To get help on what functions are available you can do: 

```sh
docker run somisana_opendrift -h
```

Or to get help for a particular function (e.g. `run_model`) you can do:

```sh
docker run somisana_opendrift run_model -h
```

To actually run the model:

```sh
docker run \
        --rm \ # remove the container running the image once it's finished running the command
        -v <your-local-run-dir>:/mnt/config_dir \ # you need to mount the local directory where you want to run the model into the container
        -v <your-forcing-dir>:/mnt/forcing \ # maybe you have a directory with your forcing files? In which case your config.py file would need to point to e.g. /mnt/forcing/your-file.nc
          somisana_opendrift \
          run_model \
            --model_type oil \
            --config_dir /mnt/config_dir \ # this reflects where you mounted <your-local-run-dir> into the container
``` 

So, the docker image allows you to run the model, do some processing, and even some plots/animations, without having to even set up a local python environment. It is intended that this approach will be useful in running the model operationally.

# Running stochastic simulations

We have added some tools for generating an ensemble of stochastic OpenDrift simulations (see `opendrift_tools/stochastic.py`). This is useful for planning purposes, for example in comparing the likely transport pathways from different potential release locations of some or other pollutant. The goal is to generate many (hundreds) of simulations, which differ only their starting time, thereby sampling the temporal variability of the environmental forcing. These simulations also only need a valid `config.py` file, set up in the same way as you would to run a single simulation. It is however useful to use wildcards (`*`) in your forcing file names in case your simulations extend beyond the temporal range of a single forcing file.

In the event that you're using the `somisana_opendrift` docker image created as above, here would be an example of how one would run a set of 50 simulations, starting on 20 Dec 2017, with the starting time of each simulation incrementing by 3.65 days:

```sh
docker run \
        --rm \
        -v <your-local-run-dir>:/mnt/run_dir \
        -v <your-forcing-dir>:/mnt/forcing \
          somisana_opendrift \
          run_stochastic \
            --model_type oil \
            --run_dir /mnt/run_dir \
            --date_start 20171220_00 \
            --run_id 1 \
            --increment_days 3.65 \
            --run_id_end 50
``` 

There are also cli functions for gridding all these runs (`grid_stochastic`), computing statistics from the ensemble of gridded outputs (`gridded_stats`) and doing a stochastic mass balance from oil spill simulations (`stochastic_massbal`). Have a look at the help for these functions on how they would be called.

# Running the model via a github action

TODO
