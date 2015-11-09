#!/bin/sh
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

emails="gmail shyfem shyfem_aux shyfem_user"

#------------------------------------------------------------------

YesNo()
{
  echo -n "$1 ? (y/n) : " | cat >&2
  read yesno
  echo "$yesno"
}

#------------------------------------------------------------------

mail="YES"
if [ "$1" = "-no_mail" ]; then
  mail="NO"
  shift
fi

file=$1

if [ $# -eq 0 ]; then
  echo "Usage: mail_shyfem.sh [-no_mail] tar-file"
  exit 1
elif [ ! -f "$file" ]; then
  echo "*** no such file: $file ...aborting"
  exit 3
fi

#------------------------------------------------------------------

echo ""								 > $tmpfile
echo "Dear All,"						>> $tmpfile
echo ""								>> $tmpfile
echo "a new shyfem file $file is available for download."	>> $tmpfile
echo "Please use the following link to download the file:"	>> $tmpfile
echo "$link"							>> $tmpfile
echo "Alternatively you can get the code directly from:"	>> $tmpfile
echo "$gitlink"							>> $tmpfile
echo ""								>> $tmpfile
echo "Release notes:"						>> $tmpfile
$fembin/extract_release.pl $FEMDIR/RELEASE_NOTES		>> $tmpfile
echo "Other information can be found in RELEASE_NOTES."		>> $tmpfile
echo ""								>> $tmpfile
echo "Best regards, Georg"					>> $tmpfile
echo ""								>> $tmpfile

#------------------------------------------------------------------

echo "Email message:"
cat $tmpfile
if [ $mail = "YES" ]; then
  answer=`YesNo "Do you want to upload and email?"`
else
  answer=`YesNo "Do you want to upload?"`
fi
[ "$answer" = "y" ] || exit 0

echo "uploading and emailing..."

#------------------------------------------------------------------

echo "uploading file $file to google drive..."
drive upload --file $file --parent $shyfemdir
status=$?
[ $status -ne 0 ] && echo "*** error uploading file" && exit 1

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

