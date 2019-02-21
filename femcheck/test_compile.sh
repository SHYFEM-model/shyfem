#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# compiles with different available compilers
#
#--------------------------------------------------------

compilers="GNU_GFORTRAN INTEL"
#compilers="GNU_GFORTRAN"

rules_arc_dir=./arc/rules
rules_dist_dir=./femcheck/rules

rules_save=$rules_arc_dir/Rules.save
rules_dist=$rules_dist_dir/Rules.dist

femdir=$( pwd )
export FEMDIR=$femdir

#--------------------------------------------------------

#trap Clean_up SIGHUP SIGINT SIGTERM
trap Clean_up 1 2 15

Clean_up()
{
  echo "Cleaning up before exiting..."
  Clean_after
  exit
}

Clean_before()
{
  rm -f *.out *.tmp
  mkdir -p $rules_arc_dir
  mv --backup=numbered ./Rules.make $rules_save
  cp $rules_dist ./Rules.make
  [ -f allstdout.txt ] && rm allstdout.txt
  [ -f allstderr.txt ] && rm allstderr.txt
}

Clean_after()
{
  rm -f *.tmp
  rm -f stdout.out stderr.out
  #cp ./Rules.make ./Rules.last		#save last Rules.make for inspection
  mv -f $rules_save ./Rules.make
  [ -f allstdout.txt ] && mv allstdout.txt allstdout.tmp
  [ -f allstderr.txt ] && mv allstderr.txt allstderr.tmp
}

SetUp()
{
  mkdir -p $rules_arc_dir $rules_dist_dir
  if [ $? -ne 0 ]; then
    echo "Cannot create directory arc... aborting"
    exit 1
  fi
  [ -f $rules_dist ] || cp ./Rules.make $rules_dist
}

WrapUp()
{
  lines=`cat allstderr.tmp | grep -v 'setting macros' | wc -l`
  echo "================================="
  echo "Final message on all compilations: "
  if [ $lines -ne 0 ]; then
    echo "  *** some errors occured in compilation..."
    echo "  (see allstderr.tmp for more details)"
  else
    echo "  no compilation errors"
  fi
  echo "================================="
}

#--------------------------------------------------------

Comp()
{
  echo "-----------------"
  Rules "$1"

  echo "start compiling in" `pwd`
  rm -f stdout.out stderr.out
  touch stdout.out stderr.out
  make cleanall > tmp.tmp 2> tmp.tmp
  make fem > stdout.out 2> stderr.tmp

  [ -f stderr.tmp ] && cat stderr.tmp | grep -v "ar: creating" > stderr.out

  lines=`cat stderr.out | wc -l`
  if [ $lines -ne 0 ]; then
    echo "*** compilation errors..."
    cat stderr.out
  else
    echo "no compilation errors"
  fi
  echo "-----------------"

  cat stdout.out >> allstdout.txt
  cat stderr.out >> allstderr.txt

  Regress
}

RulesReset()
{
  echo "resetting Rules.make file..."
  cp $rules_dist ./Rules.make
}

Rules()
{
  fembin/subst_make.pl -quiet -first "$1" Rules.make > tmp.tmp
  mv tmp.tmp Rules.make
  #fembin/subst_make.pl   "$1" Rules.make > tmp.tmp

  echo "setting macros: $1"
  echo "setting macros: $1" >> allstdout.txt
  echo "setting macros: $1" >> allstderr.txt
}

Regress()
{
  [ "$regress" = "NO" ] && return

  echo "running regression test..."
  make regress >> allstdout.txt
  cd femregress
  make status
  cd ..
  echo "finished running regression test..."
}

#--------------------------------------------------------------------

regress="NO"
if [ "$1" = "-regress" ]; then
  regress="YES"
fi

SetUp
Clean_before

for comp in $compilers
do

  echo "================================="
  echo "compiling with $comp"
  echo "================================="
  RulesReset
  Rules "FORTRAN_COMPILER=$comp"

  make compiler_version > /dev/null 2>&1

  if [ $? -ne 0 ]; then
    echo "*** compiler $comp is not available..."
    continue
  else
    echo "compiler $comp is available..."
  fi

  Comp "ECOLOGICAL=NONE GOTM=true NETCDF=false SOLVER=SPARSKIT PARALLEL_OMP=false PARALLEL_MPI=NONE"
  #Comp "ECOLOGICAL=EUTRO GOTM=false SOLVER=PARDISO"
  Comp "ECOLOGICAL=EUTRO GOTM=false"
  #Comp "ECOLOGICAL=ERSEM GOTM=true NETCDF=true SOLVER=GAUSS"
  Comp "ECOLOGICAL=NONE GOTM=true NETCDF=true SOLVER=SPARSKIT"
  Comp "ECOLOGICAL=AQUABC NETCDF=false PARALLEL_OMP=true"

  [ "$regress" = "NO" ] && continue

  Rules "ECOLOGICAL=NONE GOTM=true NETCDF=false SOLVER=SPARSKIT PARALLEL_OMP=false PARALLEL_MPI=NONE"

  Comp "COMPILER_PROFILE=SPEED PARALLEL_OMP=true"
  Comp "COMPILER_PROFILE=CHECK PARALLEL_OMP=false"
done

Clean_after
WrapUp

#--------------------------------------------------------------------

