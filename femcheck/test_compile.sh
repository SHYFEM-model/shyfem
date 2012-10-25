#!/bin/sh
#
# compiles with differetn available compilers
#
#--------------------------------------------------------

compilers="GNU_GFORTRAN INTEL"
#compilers="GNU_GFORTRAN"

#--------------------------------------------------------

#trap Clean_up SIGHUP SIGINT SIGTERM
trap Clean_up 1 2 15

Clean_up() {

  echo "Cleaning up before exiting..."
  Clean_after
  exit
}

Clean_before() {

  rm -f *.out *.tmp
  mv --backup=numbered ./Rules.make ./rules/Rules.save
  cp rules/Rules.dist ./Rules.make
  [ -f allstdout.txt ] && rm allstdout.txt
  [ -f allstderr.txt ] && rm allstderr.txt
}

Clean_after() {

  rm -f *.tmp
  rm -f stdout.out stderr.out
  mv -f ./rules/Rules.save ./Rules.make
  [ -f allstdout.txt ] && mv allstdout.txt allstdout.tmp
  [ -f allstderr.txt ] && mv allstderr.txt allstderr.tmp
}

#--------------------------------------------------------

Comp()
{
  Rules "$1"

  echo "start compiling..."
  rm -f stdout.out stderr.out
  touch stdout.out stderr.out
  make cleanall > tmp.tmp 2> tmp.tmp
  make fem > stdout.out 2> stderr.tmp

  [ -f stderr.tmp ] && cat stderr.tmp | grep -v "ar: creating" > stderr.out

  lines=`cat stderr.out | wc -l`
  if [ $lines -ne 0 ]; then
    echo "*** errors compiling..."
    cat stderr.out
  else
    echo "no compiling errors"
  fi

  cat stdout.out >> allstdout.txt
  cat stderr.out >> allstderr.txt
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

#--------------------------------------------------------------------

Clean_before

for comp in $compilers
do
  echo "================================="
  echo "compiling with $comp"
  echo "================================="
  Rules "COMPILER=$comp"

  Comp "ECOLOGICAL=NONE GOTM=true NETCDF=false SOLVER=SPARSKIT PARALLEL=false"
  #Comp "ECOLOGICAL=EUTRO GOTM=false SOLVER=PARDISO"
  Comp "ECOLOGICAL=EUTRO GOTM=false"
  Comp "ECOLOGICAL=ERSEM GOTM=true NETCDF=true SOLVER=GAUSS"
  Comp "ECOLOGICAL=AQUABC NETCDF=false PARALLEL=true"
done

Clean_after

