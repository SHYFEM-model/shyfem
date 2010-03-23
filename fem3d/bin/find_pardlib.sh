#!/bin/sh

mkldir="/opt/intel/mkl"

# Find machine type
machine=`uname -m`
if [ $machine = 'x86_64' ]; then
	machdir='em64t'
	lib1='libmkl_solver_lp64.a'
	lib2='libmkl_intel_lp64.a'
elif [ $machine = 'i686' ]; then
	machdir='32'
	lib1='libmkl_solver.a'
	lib2='libmkl_intel.a'
else
	echo "Unknown machine: $machine"
	exit 1
fi

# Libraries with the same name
lib3='libmkl_intel_thread.a'
lib4='libmkl_core.a'
lib5='libguide.a'

# Find the last version directory
versdir=`ls -rt $mkldir | tail -1`

# Name of the libraries dir
DIRLIB_MKL=$mkldir/$versdir/lib/$machdir

# String
LIBG_MKL="-L$DIRLIB_MKL $DIRLIB_MKL/$lib1 $DIRLIB_MKL/$lib2 -Wl,--start-group $DIRLIB_MKL/$lib3 $DIRLIB_MKL/$lib4 -Wl,--end-group $DIRLIB_MKL/$lib5 -lpthread"
echo "$LIBG_MKL"

exit 0
