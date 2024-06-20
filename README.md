# Adding Fe2+ and Ferrohydrite to Hommel2019a

I’ve been trying to add a new component Fe, and a new mineral, also a new reaction Ferrohydrite Fe(OH)2 = Fe + 2OH to the original system written in Hommel2019a. But I couldn’t get it to run.


To add them, I changed the code in the following files:

•	appl/icp/eicp/leocolumnproblem.hh, leo.cc, leo_column1.input

defined initial condition, boundary condition (injection) for Fe and Ferrohydrite

•	dumux/material/fluidsystem/leomin.hh 

added Fe to the primary component by setting numComponents = 9

defined the property of Fe

•	dumux/material/solidsystem/leominsolids.hh

added Ferrohydrite

defined the property of Ferrohydrite

•	dumux/material/chemistry/leocarbonicacid.hh 

added the function for calculating reaction rate for Ferrohydrite, including frdiss and frprec

added Fe to charge balance, total molarity, et al

•	dumux/material/components

added iron2.hh, ferrohydrite.hh


To compile and run columnleobox, I follow the steps below:

./dune-common/bin/dunecontrol bexec rm -r CMakeFiles CMakeCache.txt

./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all



cd ./build-cmake/appl/icp/eicp

make columnleobox



./columnleobox leo_column1.input



The compile log file sees compile.log.

The error file sees error.log.
