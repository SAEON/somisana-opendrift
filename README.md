# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used to run [OpenDrift](https://opendrift.github.io/) simulations on top of our [operational ocean modelling system](https://github.com/SAEON/somisana-croco). The [OpenDrift](https://opendrift.github.io/) documentation is very good, and one should start there if new to OpenDrift. The primary purpose of this repo is to provide tome tools to run OpenDrift simulations operationally on a server set up for this purposes e.g. in response to search and rescue or oil spill emergencies. It is also possible to clone this repo locally, in which case it provides a way of running your own local simulations with minimal setup.  

Directories in the repository: 
- `opendrift_tools`:   python functions for some preprocessing, running, postprocessing and plotting OpenDrift simulations
- `configs`:           configuration files used for running OpenDrift simulations - these files are used to define the variables used in various pre and postprocessing functions
- `.github/workflows`: github workflows for running OpenDrift operationally on a server set up for this purpose on MIMS 

# Setting up the python environment for local use and/or development

To run OpenDrift simulations locally, you need to first have the `opendrift` mamba/conda environment set up, according to the [OpenDrift install instructions](https://opendrift.github.io/install.html).

Once you've set up your opendrift environment, you can clone this repo to your local machine (do this wherever you would like the code):

`git clone git@github.com:SAEON/somisana-opendrift.git`

Then navigate to the root directory of the repo, and install the code into your `opendrift` environment:

```sh
cd somisana-opendrift
conda activate opendrift
pip install --no-deps -e .
```

# Example of running a simulation locally

To run a local simulation from this repo, all that is needed is to set up a relevant `config_*.py` file which are in the `configs` directory. These default configuration files are set up to run an operational simulation. So to run locally, you would need to copy the relevant default file to wherever you want to run the model, and then edit it accordingly - particularly the paths to the forcing files which will obviously not be relevant for your local run. For example, to run an oil spill simulation:

```sh
cp somisana-opendrift/configs/config_oil.py <your-config-dir>/config.py
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

It is often desirable to go from the Lagrangian particle output to a Eulerian field. So far this repo contains functions for calculating the surface oil thickness, stranded oil concentration, or the particle density (which would be applicable to all types of runs, not just oil spills). To call the function from your own python script:

```sh
from opendrift_tools.postprocess import grid_particles
grid_particles(...) # check out the documentation in the function for the expected inputs
```

Or again, you can do this from the command line interface (cli.py), in which case the inputs to the function are read from the config.py file:

```sh
python /home/gfearon/code/somisana-opendrift/cli.py grid_particles --config_dir <your-config-dir>
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

But if you want to personalise your plot more, it's probably best to do your own thing, perhaps based on the `plot_particles()` and `plot_gridded()` functions in `opendrift_tools/plotting.py`... 

# Running the docker image

TODO

# Running the model via a github action

TODO
