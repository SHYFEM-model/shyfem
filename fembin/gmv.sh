#!/bin/sh
#
# returns major version of gfortran
#
#-----------------------------------------------------------

prog=$( which gfortran )

if [ -n "$prog" ]; then
  $prog --version			\
		| head -1		\
                | sed -e 's/.*) *//' 	\
		| sed -e 's/\..*//'
fi

