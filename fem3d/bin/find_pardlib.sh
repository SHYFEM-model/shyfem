#!/bin/bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# Checks include directory
if [ -d $MKLROOT/include ]; then
   MKLINCLUDE=$MKLROOT/include
else
   echo "find_pardlib.sh: $MKLINCLUDE does not exist."
   exit 1
fi

# Find machine type
machine=`uname -m`
if [ $machine = 'x86_64' ]; then
        if [ -d $MKLROOT/lib/em64t ]; then
           DIRLIB_MKL=$MKLROOT/lib/em64t
           INTEL_INCL=$MKLROOT/../compiler/include/em64t
        elif [ -d $MKLROOT/lib/intel64 ]; then
           DIRLIB_MKL=$MKLROOT/lib/intel64
           INTEL_INCL=$MKLROOT/../compiler/include/intel64
        else
           echo "find_pardlib.sh: $DIRLIB_MKL does not exist."
           exit 1
        fi
	baselib='mkl_intel_lp64'

elif [ $machine = 'i686' ]; then
        if [ -d $MKLROOT/lib/32 ]; then
           DIRLIB_MKL=$MKLROOT/lib/32
           INTEL_INCL=$MKLROOT/../compiler/include/32	#to check
        elif [ -d $MKLROOT/lib/ia32 ]; then
           DIRLIB_MKL=$MKLROOT/lib/ia32
           INTEL_INCL=$MKLROOT/../compiler/include/ia32	#to check
        else
           echo "find_pardlib.sh: $DIRLIB_MKL does not exist."
           exit 1
        fi
	baselib='mkl_intel'

else
	echo "Unknown machine: $machine"
	exit 1
fi

# This is not necessary if '!$' is used in pardiso_solve
# before 'use omp_lib'
#if [ -e $INTEL_INCL/omp_lib.mod ]; then
#   cp $INTEL_INCL/omp_lib.mod ../femlib/mod
#else
#   echo "Module file omp_lib.mod not found"
#   exit 1
#fi

######
# Static linking. Warning: not sure if possible with GNU licence
#LIBG_MKL="-L$DIRLIB_MKL -I$MKLINCLUDE -Wl,--start-group $DIRLIB_MKL/lib${baselib}.a $DIRLIB_MKL/libmkl_intel_thread.a $DIRLIB_MKL/libmkl_core.a -Wl,--end-group -liomp5 -lpthread"

# Dynamic linking
LIBG_MKL="-L$DIRLIB_MKL -I$MKLINCLUDE -l${baselib} -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"

echo "$LIBG_MKL"

exit 0
