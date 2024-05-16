# freezing_branch
This MATLAB code simulates freeze induced cell dehydration, pressure and diameter changes in a tree stem.
Here is the associated publication : https://doi.org/10.1093/treephys/tpad117

The main directory contains different sub-directories that correspond to different simulations. \
Each sub-directory contains three .m files: main.m, parameters.m and dyfun.m (see explanation below), as well as data in .mat format (e.g., pressure and temperature fields, ...) produced by the simulations. \
The main directory contains experimental data, and two .m files : comparison.m and post_process.m.

Author: Cyril Bozonnet (cyril.bozonnet@inrae.fr; github: cyrilbz) 
INRAE, PIAF, Clermont-Ferrand
        
The code structure has been inspired from an existing code 
written by Isabell Graf (Konrad) and John M. Stockie
Department of Mathematics
Simon Fraser University
         
Developped using Matlab version R2018a.

The code is distributed under the CeCILL-B free software 
icense agreement.
(https://cecill.info/licences/Licence_CeCILL-B_V1-en.html)

# Code (short) description

## parameters.m:
This file contains all the parameters of the model as well as many useful functions. Everything is wrapped within the "p" structure so it is easier to access and transfer the parameters to the different parts of the code.


## dyfun.m:
This is the code that contains the model equations. It constructs the vector dy/dt, with y the vector of unknowns to be computed, using the parameters and previous (known) values of y.

## main.m:
This is the main code that runs the simulations.
It initializes y(t0), and computes y(t1) as a function of y(t0) and dy/dt using Matlab ode15s time advancement over a given timespan and with different options. It stops the time integration at regular time intervals to save the data (.mat files) so one can analyze the results on the fly. Note that it only saves the model variables and not the state variables.  

## post_process.m:
It is used to recompute the state variables, do some post-processing (e.g., diameter changes computation) and plot results. It saves post processed variables and state variables.

## comparison.m:
It is used to compare different model results (contained in different folders). Post_process must have been run on each case before comparison to be done. 

# Recommandations for users and developpers

- When modifying a state equation in dyfun.m, one must also copy the modification in post_process.m for consistency
- If one increases the solver tolerance in main.m, it might lead to convergence issues.
- Always create a new sub-directory for each new simulation: the .mat data files are overwritten at each run. The path between the post_process.m and comparison files with the sub-directories can be changed (no need to keep the same architecture as long as the main.m, parameters.m and dyfun.m files are in the same directory).
