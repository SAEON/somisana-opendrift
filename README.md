# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used to run [OpenDrift](https://opendrift.github.io/) simulations on top of our [operational ocean modelling system](https://github.com/SAEON/somisana-croco). The [OpenDrift](https://opendrift.github.io/) documentation is very good, and one should start there if new to OpenDrift. The primary purpose of this repo is to provide some tools to run OpenDrift simulations operationally on a server set up for this purposes e.g. in response to search and rescue or oil spill emergencies. It is also possible to clone this repo locally, in which case it provides a way of running your own local simulations with minimal setup.  

Directories in the repository: 
- `opendrift_tools`:   python functions for some preprocessing, running, postprocessing and plotting OpenDrift simulations
- `configs`:           configuration files used for running OpenDrift simulations - these files are used to define the variables used in various pre and postprocessing functions
- `.github/workflows`: github workflows for running OpenDrift operationally on a server set up for this purpose on DFFE's Marine Information Management System (MIMS) 

Please refer to the [wiki](https://github.com/SAEON/somisana-opendrift/wiki) for more information.

