# l1l2inv2

# DESCRIPTION
This program performs magnetic inversion with L1-L2 norm combined regularization using CDA,
via coordinate descent algorithm.

minimize (1/2) * || y - X * beta ||^2 + (1/2) * lambda * (1 - alpha) * ||beta||^2+lambda * alpha * |beta|,

where

lambda: regularization parameter

alpha: hyper-parameter which controls the strength of L1 and L2 norm regularization

beta: subsurface magnetic model

y: observed magnetic anomaly

X: transfer matrix

The optimal lambda is selected using the L-curve criterion, with a priori specification of alpha.
For more detail of the inversion, please see

https://earth-planets-space.springeropen.com/articles/10.1186/s40623-019-1052-4

# DEPENDENCIES
This program uses the following external libraries:
 * BLAS
 * GSL(Gnu Scientific Libraries)
 * openMP

Before compile this program, please install the above libraries in your system.

# INSTALLATION
Before the compilation, copy template file make.config.templ as make.config:

$ cp make.config.templ make.config

and modify the following entries in this file according to your system:

 * BLAS_LIB: BLAS library path and include dirs
 * GSL_LIB:  GSL library path and include dirs
 * OPENMP_FLG: openMP flag (e.g. -fopenmp for gcc, etc.)
 * CC: C compiler (gcc, icc, etc.)

After modifying, run make, and make install.
If you want to install this program in a specific directory, run

$ make install DESTDIR=\<root\>


# INVERSION

After finishing installation correctly, a script calc.sh should be exists in \<root\>/bin .
To perform inversion, please run this script.
You can see the usage of this script by run calc.sh with -h option:

$ \<root\>/bin/calc.sh -h

### prepare setting-file, data file and terrain data file
#### input.data
The script calc.sh searches a file named "input.data" which involves observed magnetic anomaly.
This file "input.data" has to involve site location and observed magnetic anomaly.
The requierd data format is

<EW location of obs. point(km)> <NS location(km)> <altitude(km)> <observed anomaly(nT)>

If this file "input.data" is not exists, calc.sh will abort.
Please show a sample input data file in demo/samples 

#### settings
The script calc.sh also search a file named "settings" which involves some parameter settings
such as magnetic inclination, declination, number of grid sells, etc.
For more detail, please refer a sample settings file in demo/samples/settings
(the easiest way is to copy this file to the current directory and edit it).

### terrain.data
If a file named "terrain.data", which involves gridded terrain data, exists on your current directory,
calc.sh automatically read this file and import terrain data into the inversion program.
The requierd data format of "terrain.data" is

<EW location(km)> <NS location(km)> <altitude(km)>

The dimensions and number of the terrain grid must consistent with the dimensions
and horizontal number of the grid cells that subdivide the subsurface space.

### For quite large problem

If size of the matrix X is too large, memory allocation will fail and calc.sh will terminate abnormally.
This is because calc.sh tries to store the matrix X in memory at one time.
If you treat a large X, that is, when you use very fine grid and/or large observation data,
the script calc_xmat.sh is available.
This script does not store X in memory, but in files, and read them when needed.
So, calc_xmat.sh uses only small memories, but uses large space of storage.
Before to run calc_xmat.sh, please confirm you have enough free space in your HDD.

### notice
By using calc_xmat.sh, the performance of the inversion is reduced
because of the overhead of the accessing to the storage to read the matrix X.

The scripts calc.sh and calc_xmat.sh call the following programs internally:

 * bin/l1l2inv, bin/l1l2inv_xmat: main inversion programs
 * bin/lcurve_interp: derives interpolated L-curve with the aid of the smooth b-spline
 * bin/optimal_lambda: chooses an optimal lambda using interpolated L-curve

## Output

calc.sh, and calc_xmat.sh output the models (beta) for a decreasing sequence of lambda
to an ascii file beta_path.data (or beta_path1.data).
To extract a model of specific number of iteration, use a support program bin/extract.
To see the usage of this program, please run

$ \<root\>/bin/extract -h
