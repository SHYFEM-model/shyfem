#!/bin/sh
#
# compiles with differetn available compilers
#
#--------------------------------------------------------

compilers="GNU_GFORTRAN INTEL"
#compilers="GNU_GFORTRAN"

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
    echo "errors compiling..."
    cat stderr.out
  else
    echo "no compiling errors"
  fi

  cat stdout.out >> allstdout.out
  cat stderr.out >> allstderr.out
}

Rules()
{
  fembin/subst_make.pl -quiet -first "$1" Rules.make > tmp.tmp
  mv tmp.tmp Rules.make
  #fembin/subst_make.pl   "$1" Rules.make > tmp.tmp

  echo "setting macros: $1"
  echo "setting macros: $1" >> allstdout.out
  echo "setting macros: $1" >> allstderr.out
}

#--------------------------------------------------------------------

rm -f *.out *.tmp
mv --backup=numbered ./Rules.make ./Rules.save
cp rules/Rules.test ./Rules.make

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

rm -f *.tmp
rm -f stdout.out stderr.out
mv -f ./Rules.save ./Rules.make
mv allstdout.out allstdout.tmp
mv allstderr.out allstderr.tmp

