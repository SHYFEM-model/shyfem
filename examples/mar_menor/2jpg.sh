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
# converts ps files to jpg
#
#-----------------------------------------------

ConvertJpg()
{
  file=$1
  keep=$2	#we keep that frame

  name=$( basename $file .ps )
  echo "treating: $name   keeping: $keep"

  gps -split -eps -jpg -dpi 200 $file

  CleanFiles $name $keep
}

CleanFiles()
{
  name=$1
  keep=$2

  rm -f $name.[0-9].ps
  rm -f $name.[0-9].eps

  [ -z "$keep" ] && return

  for n in $name.[0-9].jpg
  do
    if [ $n = $name.$keep.jpg ]; then
      echo "keeping $n"
    else
      echo "deleting $n"
      rm -f $n
    fi
  done
}

#-----------------------------------------------

ConvertJpg mm_hyd_31.salt.apnbath.ps 7
ConvertJpg mm_hyd_31.salt.ps
ConvertJpg mm_hyd_32.salt.ps

ConvertJpg mm_hyd_41.salt.1.apnbath.ps 7
ConvertJpg mm_hyd_41.salt.5.apnbath.ps 7
ConvertJpg mm_hyd_41.vel.1.apnbath.ps 7
ConvertJpg mm_hyd_41.vel.5.apnbath.ps 7

ConvertJpg mm_hyd_42.vel.1.apnbath.ps 7
ConvertJpg mm_hyd_42.vel.5.apnbath.ps 7
ConvertJpg mm_hyd_43.vel.1.apnbath.ps 6
ConvertJpg mm_hyd_43.vel.5.apnbath.ps 6

#-----------------------------------------------

