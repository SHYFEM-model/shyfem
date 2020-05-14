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
# uploads files to shyfem directory in google drive
#
#------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

fembin=$FEMDIR/fembin

shyfemdir="0B742mznAzyDPbGF2em5NMjZYdHc"
link="https://drive.google.com/folderview?id=$shyfemdir&usp=sharing"
gitlink="https://github.com/SHYFEM-model/shyfem"
tmpfile=tmp.tmp
#fembin=./fembin

emails="gmail shyfem_g shyfem_d shyfem_u"

#------------------------------------------------------------------

YesNo()
{
  echo -n "$1 ? (y/n) : " | cat >&2
  read yesno
  echo "$yesno"
}

Help()
{
  echo "Usage: copy_shyfem.sh [-h|-help] file"
  exit 1
}

Gversion()
{
  gver=$( gdrive -v 2> /dev/null )
  if [ -z "$gver" ]; then
    aux=$( gdrive version 2> /dev/null )
    #echo "auxxxxxxxx: $aux ----"
    gver=$( echo "$aux" | head -1 )
  fi
  #echo "status: $?"
  echo "gdrive version: $gver"
}

Clean()
{
  [ -f tmp.tmp ] && rm -f tmp.tmp
}

#------------------------------------------------------------------

file=$1

if [ $# -eq 0 ]; then
  echo "no file given..."
  Help
fi

#------------------------------------------------------------------

echo "uploading file $file to google drive $shyfemdir"
Gversion
if [ "$gver" = "gdrive v1.9.0" ]; then			#for 1.9.0
  gdrive upload --file $file --parent $shyfemdir
elif [ "$gver" = "gdrive v2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $shyfemdir $file
elif [ "$gver" = "gdrive: 2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $shyfemdir $file
else
  echo "unknown version of gdrive: $gver"
  exit 1
fi
status=$?
[ $status -ne 0 ] && echo "*** error uploading file" && exit 1

#------------------------------------------------------------------

