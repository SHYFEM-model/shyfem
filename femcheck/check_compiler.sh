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
# checks if all compilers are available
#
#------------------------------------------------------------------------

#F77
#F95
#cc

debug=YES
debug=NO

TestCompiler()
{
  [ $debug = YES ] && echo "test compiler $1"
  if [ -z "$1" ]; then
    echo "*** no compiler given: $1"
    exit 3
  fi
  command $1 > /dev/null 2>&1
  status=$?
  if [ $status -ne 0 ]; then
    echo "*** error executing compiler $1"
    exit 1
  else
    :
    [ $debug = YES ] && echo "  ...ok compiler $1"
  fi
}

for compiler
do
  TestCompiler "$compiler"
done

exit 0

