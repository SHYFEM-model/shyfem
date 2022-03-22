#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
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

subject="new SHYFEM release"
shyfemdir="0B742mznAzyDPbGF2em5NMjZYdHc"
link="https://drive.google.com/folderview?id=$shyfemdir&usp=sharing"
gitlink="https://github.com/SHYFEM-model/shyfem"
tmpfile=tmp.tmp
addresses=$FEMDIR/femcheck/emails/to_mail.txt

#------------------------------------------------------------------

YesNo()
{
  echo -n "$1 ? (y/n) : " | cat >&2
  read yesno
  echo "$yesno"
}

Help()
{
  echo "Usage: mail_shyfem.sh [-h|-help] [-no_mail] [-no_upload] tar-file(s)"
  exit 0
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

MakeLetter()
{
  echo ""							 > $tmpfile
  echo "Dear All,"						>> $tmpfile
  echo ""							>> $tmpfile
  echo "a new shyfem release is available for download."	>> $tmpfile
  echo "Please go to Github and download the new release:"	>> $tmpfile
  echo "$gitlink"						>> $tmpfile
  echo "Click on \"releases\" and choose the desired version"	>> $tmpfile
  echo ""							>> $tmpfile
  echo "Release notes:"						>> $tmpfile
  $fembin/extract_release.pl $FEMDIR/RELEASE_NOTES		>> $tmpfile
  echo "Other relevant information can be found in:"		>> $tmpfile
  echo "    RELEASE_NOTES, LOG, COMMIT, VERSION, BUG"		>> $tmpfile
  echo ""							>> $tmpfile
  echo "If you do not want to obtain these kind of messages"	>> $tmpfile
  echo "please let me know and I will delete you from this list.">> $tmpfile
  echo ""							>> $tmpfile
  echo "Best regards, Georg"					>> $tmpfile
  echo ""							>> $tmpfile
}

UploadFiles()
{
  UploadFile $1
  shift
  UploadFile $1
}

UploadFile()
{
  file=$1
  [ -n "$file" ] || return
  if [ ! -f $file ]; then
    echo "*** no such file: $file... exiting"
    return
  fi

  echo "uploading file $file to google drive..."
  Gversion

  if [ "$gver" = "gdrive v1.9.0" ]; then		#for 1.9.0
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
}

MailMessage()
{
  local ilines=0
  local file=$addresses

  if [ ! -f $file ]; then
    echo "*** cannot find file $file"
    return
  fi

  while read -r line; do
    (( ilines++ ))
    echo "$line"
  done < $file

  echo "$ilines emails found"
  answer=`YesNo "Do you want to email to these addresses?"`
  [ "$answer" = "y" ] || return

  while read -r line; do
    [[ $line = Georg* ]] || continue
    email=$( echo $line | sed -e 's/.*</</' )
    echo "sending mail to $email"
    gmutt -auto -s "$subject" -i $tmpfile $email
  done < $file
}

MailMessage0()
{
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
}

#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------

mail="YES"
upload="NO"

while [ -n "$1" ]
do
   case $1 in
        -no_mail)       mail="NO";;
        -no_upload)     upload="NO";;
        -h|-help)       Help; exit 0;;
        -*)             echo "No such option $1"; exit 1;;
        *)              break;;
   esac
   shift
done

file1=$1
file2=$2

[ $# -eq 0 ] && upload="NO"

#------------------------------------------------------------------

MakeLetter

echo "Email message:"
cat $tmpfile
if [ $upload = "YES" ]; then
  echo "Files to be uploaded:"
  echo "  $file1 $file2"
fi

echo "uploading and emailing..."

#------------------------------------------------------------------

#if [ $upload = "YES" ]; then
#  answer=`YesNo "Do you want to upload?"`
#  [ "$answer" = "y" ] && UploadFiles $file1 $file2
#fi

#if [ $mail = "YES" ]; then
#  answer=`YesNo "Do you want to email?"`
#  [ "$answer" = "y" ] && MailMessage
#fi

MailMessage

Clean

#------------------------------------------------------------------

