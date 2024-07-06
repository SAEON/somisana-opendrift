# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used to run [OpenDrift](https://opendrift.github.io/) simulations on top of our [operational ocean modelling system](https://github.com/SAEON/somisana-croco). 

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

# Running an example oil spill simulation locally

TODO

# Running an example oil spill simulation as a github action

TODO
