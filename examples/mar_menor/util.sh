#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------------------------------
# utility routines for running and plotting
#--------------------------------------------------------------------

Run()
{
  str=$1

  shyfem $str.str

  status=$?

  if [ $status -ne 99 ]; then	#simulation should finish with status 99
    echo "*** error in run... status:$status  expecting:99 ...aborting"
    exit 77
  fi
}

ElabSim()
{
  command=$1
  shift
  files=$*

  $command $files
  CheckStatus "$command" $?
}

CheckFiles()
{
  for file
  do
    if [ ! -f $file ]; then
      echo "file $file not existing ...aborting"
      exit 71
    fi
  done
}

CleanFiles()
{
  for file
  do
    if [ -f $file ]; then
      rm $file
    fi
  done
}

CheckStatus()
{
  prog=$1
  status=$2
  expected=$3
  [ "$expected" = "" ] && expected=0

  if [ $status -ne $expected ]; then
    echo "error running $prog  status:$status  expecting:$expected ...aborting"
    exit 73
  fi
}

#--------------------------------------------------------------------

PlotTsVel()
{
  files=$*

  gp -t "Current speed" -tx "time" -ty "current speed [m/s]" $files
  CheckStatus gp $?
  mv out.ps $sim.speed.ps
}

PlotTsZeta()
{
  files=$*

  gp -t "Water level" -tx "time" -ty "water level [m]" $files
  CheckStatus gp $?
  mv out.ps $sim.zeta.ps
}

PlotTsSalt()
{
  files=$*

  gp -t "Salinity" -tx "time" -ty "salinity [psu]" $files
  CheckStatus gp $?
  mv out.ps $sim.salt.ps
}

PlotMapVel()
{

  apn=$1
  level=$2

  if [ -n "$level" ]; then
    plevel="-layer $level"
  else
    plevel=""
    level=0
  fi

  shyplot -varid 2 -freq 32 $plevel $sim.hydro.shy $apn.str
  CheckStatus shyplot $?
  mv plot.ps $sim.vel.$level.$apn.ps
}

PlotMapSalt()
{

  apn=$1
  level=$2

  if [ -n "$level" ]; then
    plevel="-layer $level"
  else
    plevel=""
    level=0
  fi

  shyplot -varid 11 -freq 32 $plevel $sim.ts.shy $apn.str
  CheckStatus shyplot $?
  mv plot.ps $sim.salt.$level.$apn.ps
}

#--------------------------------------------------------------------

