#!/bin/sh
#
#--------------------------------------------------

YesNo()
{
  echo -n "Ok ? (y/n) : " | cat >&2
  read yesno
  echo "$yesno"
}

#--------------------------------------------------

fem=$1

if [ $# -ne 1 ]; then
  echo "Usage: pubfem.sh fem-file"
  exit 1
elif [ ! -f $fem ]; then
  echo "no such file: $fem"
  exit 1
fi

#--------------------------------------------------

who="gmail chris debora fra bajo michol cucco"
who="gmail"
letter=~/fem/publish/letter.txt

vers=`basename $fem '.tar.gz'`
subject=$vers

echo "file:    $fem"
echo "who:     $who"
echo "subject: $subject"
#echo "letter:  $letter"

ok=`YesNo`
[ "$ok" = "y" ] || exit 0

#--------------------------------------------------

echo "... copying file"

cp $fem ~
scp $fem model@150.178.42.69:SHYFEM/.

#--------------------------------------------------

echo "... sending letter"

gmutt -s "$subject" $who <<EOI

   A new version of SHYFEM has been published.

   You can find it in:

        - lagoon (georg@150.178.42.186) in directory /home/georg
        - ftp site (model@150.178.42.69) in directory SHYFEM

EOI

echo "$vers has been published"

#--------------------------------------------------
