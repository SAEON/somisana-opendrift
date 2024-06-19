# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used to run [OpenDrift](https://opendrift.github.io/) simulations on top of our [operational ocean modelling system](https://github.com/SAEON/somisana-croco). 

# Setting up the python environment for local use and/or development

To run your own OpenDrift simulations, you need to first have the `opendrift` conda environment set up, according to the [OpenDrift install instructions](https://opendrift.github.io/install.html).

Additionally, to easily access the code in this repo, you'll want to install the python code in your `opendrift` environment:

```sh
conda activate opendrift
pip install --no-deps -e .
```

# Running an example oil spill simulation locally

TODO

# Running an example oil spill simulation as a github action

TODO
