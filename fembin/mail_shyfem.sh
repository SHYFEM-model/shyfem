#!/bin/sh
#
# uploads files to shyfem directory in google drive and sends mail
#
#------------------------------------------------------------------

file=$1

if [ $# -eq 0 ]; then
  echo "Usage: mail_shyfem.sh tar-file"
  exit 1
elif [ ! -f "$file" ]; then
  echo "*** no such file: $1 ...aborting"
  exit 3
fi

subject="new SHYFEM file $file"
shyfemdir="0B742mznAzyDPbGF2em5NMjZYdHc"
link="https://drive.google.com/folderview?id=$shyfemdir&usp=sharing"
tmpfile=tmp.tmp

echo "uploading file $file to google drive..."
drive upload --file $file --parent $shyfemdir
status=$?
[ $status -ne 0 ] && echo "*** error uploading file" && exit 1

echo "sending mail..."
echo "Dear All,"						 > $tmpfile
echo "a new shyfem file $file is available for download."	>> $tmpfile
echo "Please use the following link to download the file:"	>> $tmpfile
echo "$link"							>> $tmpfile
echo "Best regards, Georg"					>> $tmpfile

#gmutt -auto -s "$subject" -i $tmpfile gmail
gmutt -auto -s "$subject" -i $tmpfile shyfem

