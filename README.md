# l1l2inv2

## DESCRIPTION
This program performs magnetic inversion with L1-L2 norm combined regularization.

## INSTLLATION
This program uses the following libraries
 * BLAS
 * GSL(Gnu Scientific Libraries)
 * openMP

Before install this program, please install the above libraries.

After that, open Makefile on the root dir, and modify the following entries according to your system:

 * BLAS_LIB: BLAS library, library-path and include dir
 * GSL_LIB:  GSL library, library-path and include dir
 * OPENMP_FLG: openMP flag

 * CC: C compiler
 * MAKE: make command
 * INSTALL: install command
 * AR: archiver

After modifying Makefile, run make, and make install.
If you want to install this programs in a specific directory, please run

$ make install DESTDIR=\<root\>


# RUN

After finishing installation correctly, a script calc.sh may be exists in \<root\>/bin .
To perform inversion, run this script.
To see usage of this script, run

$ \<root\>/bin/calc.sh

or

$ \<root\>/bin/calc.sh -h

# For quite large probrem

If size of the matrix X is too large, memory allocation will fail and calc.sh will terminate.
This is because calc.sh stores matrix X in memory at one time.
If you treat a large X, that is, when you use very fine grid and/or large observation data, the script calc_xmat.sh is available.
This script does not store X in memory, but in files, and read these files when needed.
So, calc_xmat.sh uses only small memories, but uses large space of storage (such as HDD),
and performance of inversion is reduced because of the overhead of accessing the storage to read the matrix X.
