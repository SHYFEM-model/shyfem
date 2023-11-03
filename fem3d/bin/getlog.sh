#!/bin/sh
#
# restore log for a file from an older commit
#
#------------------------------------------------------

if [ $# -ne 3 ]; then
  echo "Usage: getlog.sh commit old_file new_file"
  exit 0
fi

commit=$1
old_file=$2
new_file=$3


git co $commit -- $old_file
[ $? -ne 0 ] && exit 1

exit 0

mkdir -p tmp
cp $new_file tmp
[ $? -ne 0 ] && exit 1

#git mv $old_file $new_file

echo "   $old_file ->  $new_file"

#------------------------------------------------------

