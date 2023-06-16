# freezing_branch
This MATLAB code simulates freeze induced cell dehydration, pressure and diameter changes in a tree stem. \


The main directory contains different sub-directories that correspond to different simulations. \
Each sub-directory contains three .m files: main.m, parameters.m and dyfun.m (see explanation below), as well as data (e.g., pressure and temperature fields, ...) produced by the simulations. \
The main directory contains experimental data, and two .m files : comparison.m and post_process.m.

## parameters.m:
This file contains all the parameters of the model as well as many useful functions. Everything is wrapped within the "p" structure so it is easier to access and transfer the parameters to the different parts of the code.
