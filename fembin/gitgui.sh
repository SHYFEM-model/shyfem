#!/bin/bash
#
# GUI for "git status" to diff and edit files
#
#----------------------------------------------------------------

files=$( git s | grep modified: | sed -e 's/\s*modified:\s*//' )

n=0
for file in $files
do
  n=$(( n+1 ))
  list="$list $n $file"
  #echo $n $file
  array[$n]=$file
done

w=50
h=$(( n+10 ))
buttons="--extra-button --extra-label edit \
		--ok-label diff --cancel-label exit"

choice=0

#------------------------------------------------------------
# diff:0  exit:3  exit:1
#------------------------------------------------------------

while :
do

  default="--default-item $choice"

  choice=$( dialog $buttons $default \
		--menu "Choose file:" \
		$h $w $n $list 3>&2 2>&1 1>&3 )

  status=$?
  clear

  #echo "exit status = $status"
  #echo "choice = $choice"

  if [ $status -eq 1 ]; then
    break
  elif [ $status -eq 0 ]; then
    file=${array[$choice]}
    echo "diffing $file"
    gitdiff $file
  elif [ $status -eq 3 ]; then
    file=${array[$choice]}
    echo "editing $file"
    vi $file
  fi

done

#------------------------------------------------------------

