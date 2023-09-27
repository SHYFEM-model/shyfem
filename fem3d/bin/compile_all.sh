#!/bin/sh
#
# compiles several files of shyfem to check dependencies
#
#----------------------------------------------------

base=~/shyfem
act=$( pwd )
compile=$base/fem3d/bin/compile_test.sh
compiler=gfortran
tmpdir=tmp_compile
libdir=$act/tmp_lib
lib=$libdir/libcompile.a
libcommand=
libcommand="-L$libdir -lcompile"

#----------------------------------------------------

ExistsExtension()
{
  ls *.$1 > /dev/null 2>&1 || return 1
  #echo "extension $1 does exist"
  return 0
}

DeleteExtension()
{
  ls *.$1 > /dev/null 2>&1 || return
  rm -f *.$1
}

DeleteFile()
{
  [ -f $1 ]  && rm -f $1
}

DeleteFiles()
{
  DeleteExtension mod
  DeleteExtension o
  DeleteFile a.out
  DeleteFile $file
}

#----------------------------------------------------

#files=$( ls a*.f )
files=$( ls $* )
rm -f LIST

mkdir -p $libdir
mkdir -p $tmpdir
cd $tmpdir

touch $lib

echo "       end" > main_dummy.f
$compiler -c main_dummy.f

for file in $files
do
  name=${file%.*}
  ln -fs $act/$file 
  logfile=$name.log
  echo "compiling $file" > $logfile
  $compiler -c $file $libcommand >> $logfile 2>&1
  status=$?
  if [ $status -ne 0 ]; then
    echo "$status   $file" 
    DeleteFiles
    continue
  fi
  $compiler main_dummy.f $file $libcommand >> $logfile 2>&1
  status=$?
  echo "$status   $file  ($name)"
  if [ $status -eq 0 ]; then
    echo "$status   $file" >> $act/LIST
    #ls
    if  ExistsExtension mod; then
       mv -f *.mod $libdir
    fi
    ar rvs $lib $name.o > /dev/null
    [ -f $name.o ] && mv $name.o $libdir
  fi
  DeleteFiles
done

rm -f main_dummy.f
rm -f *.log
cd $act
rmdir $tmpdir

echo "objects and mod files are in $libdir"

#----------------------------------------------------

