#!/bin/sh
#
# links petsc mod files for usage
#
#--------------------------------------------------------------

petsc_dir=/usr/lib/petscdir/petsc3.10/x86_64-linux-gnu-real
shyfem_dir=~/shyfem
mod_dir=$shyfem_dir/femlib/mod

if [ ! -d $petsc_dir ]; then
  echo "*** cannot find directory $petsc_dir"
  exit 1
fi

Link()
{
  ln -fs $petsc_dir/include/$1 $mod_dir/$1
}

#--------------------------------------------------------------

Link petscdm.mod
Link petscdmlabel.mod
Link petscksp.mod
Link petscmat.mod
Link petscvec.mod


