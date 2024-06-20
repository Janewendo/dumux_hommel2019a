#!/bin/bash

# download Dune core modules
git clone https://gitlab.dune-project.org/core/dune-common.git
cd dune-common
git checkout releases/2.6
cd ..
git clone https://gitlab.dune-project.org/core/dune-geometry.git
cd dune-geometry
git checkout releases/2.6
cd ..
git clone https://gitlab.dune-project.org/core/dune-grid.git
cd dune-grid
git checkout releases/2.6
cd ..
git clone https://gitlab.dune-project.org/core/dune-istl.git
cd dune-istl
git checkout releases/2.6
cd ..
git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
cd dune-localfunctions
git checkout releases/2.6
cd ..


### DUMUX
git clone -b releases/3.1 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
git clone https://git.iws.uni-stuttgart.de/dumux-pub/hommel2019a.git hommel2019a

### run dunecontrol
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts all
