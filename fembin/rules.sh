#!/bin/sh
#
# changes Rules.make file according to given instructions
#
#----------------------------------------------------------

ShowTargets()
{
  echo "possible targets:"
  echo "  std		set standard options"
  echo "  netcdf	set netcdf"
  echo "  mpi		set mpi on nodes"
  echo "  omp		set omp"
  echo "  petsc		set petsc solver"
  echo "  intel		set intel compiler"
  echo "  ww3		set ww3 wave model"
  echo "local targets:"
  echo "  nemunas	use nemunas server settings"
}

SetMacro()
{
  echo "setting $1 to $2"
  cat Rules.make | sed -E "s|^$1|#$1|" > rules.aux
  cat rules.aux | sed -E "s|^#$1 *= *$2 *$|$1 = $2|" > Rules.make
}

SetNewMacro()
{
  echo "setting $1 to $2"
  cat Rules.make | sed -E "s|^$1|#$1|" > rules.aux
  cat rules.aux | sed -E "s|^#$1 *= *$|$1 = $2|" > Rules.make
}

CleanUp()
{
  [ -f rules.aux ] && rm -f rules.aux
}

#----------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: rules.sh [-list] target(s)"
  exit 1
fi

if [ $1 = "-list" ]; then
  ShowTargets
  exit 0
fi

#----------------------------------------------------------

for target in $*
do
  echo "preparing for target $target"

  if [ $target = "std" ]; then
    make rules_dist
  elif [ $target = "netcdf" ]; then
    SetMacro NETCDF true
  elif [ $target = "intel" ]; then
    SetMacro FORTRAN_COMPILER INTEL
    SetMacro C_COMPILER INTEL
  elif [ $target = "mpi" ]; then
    SetMacro PARALLEL_MPI NODE
    SetMacro PARTS METIS
    SetMacro NETCDF true
  elif [ $target = "omp" ]; then
    SetMacro PARALLEL_OMP true
  elif [ $target = "petsc" ]; then
    SetMacro PARALLEL_MPI NODE
    SetMacro PARTS METIS
    SetMacro SOLVER PETSC
    SetMacro NETCDF true
  elif [ $target = "ww3" ]; then
    SetMacro PARALLEL_MPI NODE
    SetMacro PARTS PARMETIS
    SetMacro SOLVER PETSC
    SetMacro NETCDF true
    SetMacro WW3 true
  elif [ $target = "nemunas" ]; then
    SetMacro PARALLEL_OMP true
    SetMacro NETCDF true
    SetNewMacro NETCDFDIR /opt/sw
  else
    echo "target not recognized: $target"
    exit 3
  fi

done

CleanUp

exit 0

#----------------------------------------------------------

