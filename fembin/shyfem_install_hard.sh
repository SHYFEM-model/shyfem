#!/bin/sh
#
# installs shyfem through hard references
#
#-------------------------------------------------

if [ ! -f VERSION ]; then
  echo "this routine must be run in main fem directory... aborting"
  exit 1
fi

if [ "$1" = "-reset" ]; then
  option="-reset"
  extension="ggu_reset"
  echo "resetting: $option"
else
  option=""
  extension="ggu_hard"
  echo "installing: $option"
fi

dir=`pwd`

echo "using directory: $dir"

cd fembin

files=`ls -d *`

for file in $files
do
  shyfem_install_hard.pl $option $dir $file > tmp.tmp
  status=$?
  echo "status: $status  $file"
  if [ $status -ne 0 ]; then
    new=$file.$extension	# for test
    new=$file			# for real
    echo "changed... creating new file version: $file $new"
    touch -r $file tmp.tmp
    chmod +x tmp.tmp
    mv -f tmp.tmp $new
  else
    rm -f tmp.tmp
  fi
done

