#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# uploads files to shyfem directory in google drive and sends mail
#
#------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

fembin=$FEMDIR/fembin

subject="new SHYFEM file $file"
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
  echo "Usage: mail_shyfem.sh [-h|-help] [-no_mail] [-no_upload] tar-file"
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

mail="YES"
if [ "$1" = "-no_mail" ]; then
  mail="NO"
  shift
fi

upload="YES"
if [ "$1" = "-no_upload" ]; then
  upload="NO"
  shift
fi

file1=$1
file2=$2

if [ $# -eq 0 ]; then
  if [ $upload = "YES" ]; then
    Help
  fi
elif [ $1 = '-h' -o $1 = '-help' ]; then
  Help
elif [ ! -f "$file1" ]; then
  echo "*** no such file: $file1 ...aborting"
  exit 3
elif [ -n "$file2" -a ! -f "$file2" ]; then
  echo "*** no such file: $file2 ...aborting"
  exit 3
fi

#------------------------------------------------------------------

echo ""								 > $tmpfile
echo "Dear All,"						>> $tmpfile
echo ""								>> $tmpfile
echo "a new shyfem release $file is available for download."	>> $tmpfile
echo "Please use the following link to download the file:"	>> $tmpfile
echo "$link"							>> $tmpfile
echo "Alternatively you can get the code directly from:"	>> $tmpfile
echo "$gitlink"							>> $tmpfile
echo "Click on \"releases\" and choose the desired version"	>> $tmpfile
echo ""								>> $tmpfile
echo "Release notes:"						>> $tmpfile
$fembin/extract_release.pl $FEMDIR/RELEASE_NOTES		>> $tmpfile
echo "Other relevant information can be found in:"		>> $tmpfile
echo "    RELEASE_NOTES, LOG, COMMIT, VERSION, BUG"		>> $tmpfile
echo ""								>> $tmpfile
echo "If you do not want to obtain these kind of messages"	>> $tmpfile
echo "please let me know and I will delete you from this list."	>> $tmpfile
echo ""								>> $tmpfile
echo "Best regards, Georg"					>> $tmpfile
echo ""								>> $tmpfile

#------------------------------------------------------------------

echo "Email message:"
cat $tmpfile
if [ $upload = "YES" ]; then
  echo "Files to be uploaded:"
  echo "  $file1"
  [ -n "$file2" ] && echo "  $file2"
  echo ""
fi

if [ $upload = "NO" ]; then
  answer=`YesNo "Do you want to email?"`
elif [ $mail = "NO" ]; then
  answer=`YesNo "Do you want to upload?"`
else
  answer=`YesNo "Do you want to upload and email?"`
fi
[ "$answer" = "y" ] || exit 0

echo "uploading and emailing..."

#------------------------------------------------------------------

if [ $upload = "YES" ]; then

echo "uploading file $file to google drive..."
Gversion
if [ "$gver" = "gdrive v1.9.0" ]; then			#for 1.9.0
  gdrive upload --file $file1 --parent $shyfemdir
  [ -f $file2 ] && gdrive upload --file $file2 --parent $shyfemdir
elif [ "$gver" = "gdrive v2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $shyfemdir $file1
  [ -f $file2 ] && gdrive upload  --parent $shyfemdir $file2
elif [ "$gver" = "gdrive: 2.1.0" ]; then		#for 2.1.0
  gdrive upload  --parent $shyfemdir $file1
  [ -f $file2 ] && gdrive upload  --parent $shyfemdir $file2
else
  echo "unknown version of gdrive: $gver"
  exit 1
fi
status=$?
[ $status -ne 0 ] && echo "*** error uploading file" && exit 1

fi

#------------------------------------------------------------------

[ "$mail" = "NO" ] && exit 0

for email in $emails
do
  echo ""
  echo "using email $email"
  mutt -A $email
  answer=`YesNo "Do you want to email to these addresses?"`
  [ "$answer" = "y" ] || continue
  echo "sending mail to $email..."
  gmutt -auto -s "$subject" -i $tmpfile $email
done

#------------------------------------------------------------------

Clean

#------------------------------------------------------------------

