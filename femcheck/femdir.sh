#!/bin/sh

FEMDIR=${FEMDIR:-$HOME/fem}
echo "actual femdir = $FEMDIR"

#-------------------------------------------------------------

SetDir()
{

  pushd $1      > /dev/null
  local tmp=`pwd`
  popd          > /dev/null

  echo $tmp
}

IsFemDir()
{
  local vers=$1/fem3d/VERSION
  if [ -f $vers ]; then
    head -1 $vers | sed -e 's/  */ /g' | cut -f4 -d" "
  else
    echo ""
  fi
}

#-------------------------------------------------------------

if [ $# -eq 0 ]; then
  :
elif [ $1 = "-r" ]; then
  FEMDIR=${FEMINSTALL:-$HOME/fem}
  PATH=$FEMDIR/bin:$PATH
  echo "new femdir = $FEMDIR"
else
  version=`IsFemDir $1`
  if [ -n "$version" ]; then
    FEMDIR=`SetDir $1`
    PATH=$FEMDIR/bin:$PATH
    echo "new femdir = $FEMDIR"
  else
    echo "*** Error: $1 is not a FEM directory. Cannot change."
  fi
fi

hash -r
export FEMDIR

