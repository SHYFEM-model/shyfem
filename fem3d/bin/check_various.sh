#!/bin/bash
#
#----------------------------------------------------

dirs="external shyfem tools utils"

#----------------------------------------------------

GetSubdirs()
{
  du | sed -e 's/^[0-9]* *//'
}

Check()
{
  #CheckH
  #CheckNoObjs
  CheckCopyright
}

CheckCopyright()
{
  files=$( ls *.f *.f90 Makefile README 2> /dev/null )
  for file in $files
  do
    grep  'Copyright' $file > /dev/null
    status=$?
    if [ $status -ne 0 ]; then
      echo "    $status $file"
    fi
  done
}

CheckNoObjs()
{
  ffiles=$( ls *.f 2> /dev/null )
  for ffile in $ffiles
  do
    name=$( basename $ffile .f )
    obj=$name.o
    if [ ! -f $obj ]; then
      echo "    no object for $ffile"
    fi
  done

  ffiles=$( ls *.f90 2> /dev/null )
  for ffile in $ffiles
  do
    name=$( basename $ffile .f90 )
    obj=$name.o
    if [ ! -f $obj ]; then
      echo "    no object for $ffile"
    fi
  done
}

CheckH()
{
  hfiles=$( ls *.h 2> /dev/null )
  if [ -n "$hfiles" ]; then
    echo $hfiles
  fi
}

#----------------------------------------------------

cd ~/shyfemcm/shyfemcm/src/

for dir in $dirs
do
  echo $dir
  cd $dir
  thisdir=$( pwd )
  subdirs=$( GetSubdirs )
  for subdir in $subdirs
  do
    echo "  $subdir"
    #[[ "$subdir" =~ femutil ]] && continue
    cd $subdir
    Check
    cd $thisdir
  done
  cd ..
done

#----------------------------------------------------

