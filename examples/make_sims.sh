#!/bin/sh
#
#-------------------------------------------

run()
{
  str=$1
  echo "===   running $str   ==="
  shyfem $str

  if [ $? -eq 99 ] ; then
    echo "===   shyfem $str completed successfully   ==="
  else
    echo "===   ERROR running shyfem $str   ==="
    exit $?
  fi
}

#-------------------------------------------

run_all() {
  for str in `ls mm_hyd_*.str` ; do
    run $str
  done
}

#-------------------------------------------

run_all
