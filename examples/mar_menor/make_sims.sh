#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# this runs all examples
#
#-------------------------------------------

RunShell()
{
  prog=$1
  expected=$2
  [ "$expected" = "" ] && expected=0

  ./$prog
  status=$?

  if [ $status -ne $expected ]; then
    echo "status = $status  expected = $expected"
    echo "*** error executing script $prog  ...aborting"
    exit 1
  fi
}

#-------------------------------------------

RunShell mm_hyd_00.sh
RunShell mm_hyd_01.sh
RunShell mm_hyd_02.sh
RunShell mm_hyd_03.sh

RunShell mm_hyd_11.sh
RunShell mm_hyd_12.sh
RunShell mm_hyd_13.sh

RunShell mm_hyd_21.sh
RunShell mm_hyd_22.sh 77	#here we expect an error
RunShell mm_hyd_22a.sh
RunShell mm_hyd_22b.sh
RunShell mm_hyd_23.sh

RunShell mm_hyd_31.sh
RunShell mm_hyd_32.sh

RunShell mm_hyd_41.sh
RunShell mm_hyd_42.sh
RunShell mm_hyd_43.sh

#-------------------------------------------

