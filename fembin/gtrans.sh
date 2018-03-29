#!/bin/sh
#
# uploads files to shyfem directory in google drive and sends mail
#
#------------------------------------------------------------------

shyfemdir="0B742mznAzyDPbGF2em5NMjZYdHc"
transdir="0B742mznAzyDPXy02V2h2a1Mwc2M"

#------------------------------------------------------------------

Help()
{
  echo "Usage: gtrans.sh [-h|-help] file"
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

#------------------------------------------------------------------

file=$1

if [ $# -eq 0 ]; then
  Help
elif [ $1 = '-h' -o $1 = '-help' ]; then
  Help
elif [ ! -f "$file" ]; then
  echo "*** no such file: $file1 ...aborting"
  exit 3
fi

#------------------------------------------------------------------

echo "uploading file $file to google drive..."
Gversion
if [ "$gver" = "gdrive v1.9.0" ]; then			#for 1.9.0
  gdrive upload --file $file --parent $transdir
elif [ "$gver" = "gdrive v2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $transdir $file
elif [ "$gver" = "gdrive: 2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $transdir $file
else
  echo "unknown version of gdrive: $gver"
  exit 1
fi
status=$?
[ $status -ne 0 ] && echo "*** error uploading file" && exit 1

#------------------------------------------------------------------

