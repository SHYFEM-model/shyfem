#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#---------------------------------------------------------------

shydir=$HOME/shyfem
copydir=$shydir/femcheck/copyright
actdir=$( pwd )

errors=0

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

CheckExeType()
{
  # this checks executuable files
  #
  # checks types of executable files
  # finds executable files that are no scripts
  # finds scripts that are not executable
  # checks if all executable files have copyright

  echo "================================================"
  echo "--- CheckExeType: checking executable files"
  echo "================================================"

  echo "--- printing all file types that are executable"

  files=$( findf -x '*' | grep -v '/tmp/' | grep -v '/arc/' )
  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sed -e 's/ELF 64-bit LSB executable.*/ELF 64-bit LSB executable/' \
	| sed -e 's/ELF 64-bit LSB relocatable.*/ELF 64-bit LSB relocatable/' \
	| sed -e 's/GIF image.*/GIF image/' \
	| sed -e 's/.*perl.*script.*/Perl script/' \
	| sort -u
  fi

  echo "--- printing files that are executable but not scripts"

  for file in $files
  do
    exec=$( file $file | grep "ELF 64-bit LSB" )
    [ $? -eq 0 ] && continue
    #echo $file
    first=$( head -1 $file )
    if [[ ! $first =~ '#!/'.* ]]; then
      errors=$(( errors + 1 ))
      echo "*** $file is not a script"
    fi
  done
  [ $errors -ne 0 ] && exit 1

  echo "--- printing scripts that are not executable"

  files=$( findf '*' | grep -v '/tmp/' | grep -v '/arc/' )

  for file in $files
  do
    [ -d $file ] && continue
    first=$( head -1 $file )
    if [[ $first =~ '#!/'.* ]]; then
      if [ ! -x $file ]; then
	errors=$(( errors + 1 ))
        echo "*** $file is a script but is not executable"
      fi
    fi
  done
  [ $errors -ne 0 ] && exit 1

  echo "--- printing scripts that have no copyright"

  files=$( findf '*' | grep -v '/tmp/' | grep -v '/arc/' \
		| grep -v /.git/ | grep -v /GD/ \
		| grep -v /Mail-Sender | grep -v /femersem/ )

  HandleCopyright "#" script
}

CheckTexType()
{
  # this checks tex files

  echo "================================================"
  echo "--- CheckTexType: checking tex files"
  echo "================================================"

  echo "--- printing all file types that are tex files"

  files=$( findf '*.tex' | grep -v '/tmp/' | grep -v '/arc/' )
  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright "%" tex
}

CheckStrType()
{
  # this checks tex files

  echo "================================================"
  echo "--- CheckStrType: checking str files"
  echo "================================================"

  echo "--- printing all file types that are str files"

  files=$( findf '*.str' | grep -v '/tmp/' | grep -v '/arc/' )
  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright "#" text
}

CheckCType()
{
  # this checks tex files

  echo "================================================"
  echo "--- CheckCType: checking c files"
  echo "================================================"

  echo "--- printing all file types that are c files"

  files=""
  MakeFilesFromExt c h
  FilterFiles /tmp/ /arc/ 
  FilterDirs femgotm '.*.h'
  FilterDirs femersem '.*.h'
  FilterDirs fem3d '.*.h'
  FilterDirs femplot '.*.h'

  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright " \*" c
}

CheckFortranType()
{
  # this checks tex files

  echo "================================================"
  echo "--- CheckFortranType: checking fortran files"
  echo "================================================"

  echo "--- printing all file types that are fortran files"

  files=""
  MakeFilesFromExt f F90 f90 F h
  #MakeFilesFromExt  F90 
  FilterFiles /tmp/ /arc/ 
  FilterDirs femgotm '.*.F90' '.*.h'
  FilterDirs femersem '.*.F90' '.*.f90' '.*.h'
  FilterDirs grid '.*.h'
  FilterDirs mesh '.*.h'
  FilterDirs hcbs '.*.h'
  FilterDirs post '.*.h'

  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright "!" fortran
}

CheckSpecialType()
{
  # this checks tex files

  echo "================================================"
  echo "--- CheckSpecialType: checking special files"
  echo "================================================"

  echo "--- printing all file types that are special files"

  special="Makefile makefile README Rules.make Include.make"

  files=""
  for file in $special
  do
    aux=$( findf "$file" | grep -v '/tmp/' | grep -v '/arc/' \
		| grep -v /Mail-Sender | grep -v /femersem/ )
    files="$files $aux"
  done

  #echo "files: $files"
  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright "#" text
}

CheckAllType()
{
  CheckExeType
  CheckTexType
  CheckStrType
  CheckCType
  CheckFortranType
  CheckSpecialType
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

MakeFilesFromName()
{
  files=""
  for name
  do
    aux=$( findf $name )
    files="$files $aux"
  done
}

MakeFilesFromExt()
{
  for ext
  do
    aux=$( findf '*.'$ext )
    files="$files $aux"
  done
}

FilterFiles()
{
  for pattern
  do
    aux=$( echo $files | tr ' ' '\n' | grep -v $pattern )
    files=$aux
  done
}

FilterDirs()
{
  dir=$1; shift

  thisdir=$( echo $actdir | sed -e 's/.*\///' )
  if [ $thisdir = $dir ]; then
    dir="."
  else
    dir="/$dir"
  fi

  #echo "dir: $actdir $thisdir $dir"

  for pattern
  do
    aux=$( echo $files | tr ' ' '\n' | grep -v $dir/$pattern )
    files=$aux
  done
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

PrintFileType()
{
  findf '*' \
	| grep -v '/tmp/' | grep -v '/arc/' \
	| xargs file -b \
	| sed -e 's/ELF 64-bit LSB executable.*/ELF 64-bit LSB executable/' \
	| sed -e 's/ELF 64-bit LSB relocatable.*/ELF 64-bit LSB relocatable/' \
	| sed -e 's/GIF image.*/GIF image/' \
	| sed -e 's/.*perl.*script.*/Perl script/' \
	| sort -u
}

FindFileType()
{
  findf '*' \
	| grep -v '/tmp/' | grep -v '/arc/' \
	| xargs file \
	| grep "$find_type"
}

HandleCopyright()
{
  c=$1			#comment character
  type=$2
  errors=0

  for file in $files
  do
    [ -d $file ] && continue
    if [ "$type" = "script" ]; then
      first=$( head -1 $file )
      [[ ! $first =~ '#!/'.* ]] && continue
    fi
      error=0
      #head -50 $file | grep -E "^$c\s+This file is part of SHYFEM." > /dev/null
      head -50 $file | grep -E "^..\s*This file is part of SHYFEM." > /dev/null
      [ $? -ne 0 ] && error=$(( error + 1 ))
      #head -50 $file | grep -E "^$c\s+Copyright \(C\)" > /dev/null
      head -50 $file | grep -E "^..\s*Copyright \(C\)" > /dev/null
      [ $? -ne 0 ] && error=$(( error + 10 ))
      #echo "-------------- $file $error"
      if [ $error -eq 11 ]; then
        if [ $write = "YES" ]; then
          echo "*** $file has no copyright... inserting"
	  $copydir/include_copyright.sh -type $type $file
	else
          echo "*** $file has no copyright..."
	  errors=$(( errors + 1 ))
	fi
      elif [ $error -ne 0 ]; then
        echo "*** $file has damaged copyright"
	errors=$(( errors + 1 ))
      fi
  done

  if [ $errors -gt 0 -a $write = "NO" ]; then
    echo "$errors files have no or damaged copyright... use -write to insert"
    exit 1
  fi
}

ShowCopyright()
{
  errors=0
  show_copy="YES"

  files=$( findf '*' )

  for file in $files
  do
    [ -d $file ] && continue
      error=0
      head -50 $file | grep -E "^.\s+This file is part of SHYFEM." > /dev/null
      [ $? -ne 0 ] && error=$(( error + 1 ))
      head -50 $file | grep -E "^.\s+Copyright \(C\)" > /dev/null
      [ $? -ne 0 ] && error=$(( error + 10 ))
      #echo "-------------- $file $error"
      if [ $error -eq 0 ]; then
        echo "$file has copyright"
      fi
  done
}

ShowStats()
{
  errors=0
  show_copy="YES"

  dirs=$( findf -d '*' )

  for dir in $dirs
  do
      continue
      git ls-files --error-unmatch $dir > /dev/null 2>&1
      status=$?
      if [ $status -eq 0 ]; then
        echo "$dir is in git"
      else
        echo "$dir is not in git"
      fi
  done

  for dir in $dirs
  do
   for pattern in '.*/arc$' '.*/tmp$'
   do
    #echo pattern: "$pattern"
    if [[ "$dir" =~ $pattern ]]; then
      git ls-files --error-unmatch $dir > /dev/null 2>&1
      status=$?
      if [ $status -eq 0 ]; then
        echo "*** $dir is in git"
      else
	:
        echo "$dir is not in git"
      fi
    fi
   done
  done
      git ls-files info

  files=$( findf '*' )

  $copydir/copy_stats.pl $files
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

ErrorOption()
{
  echo "*** no such option: $1"
}

Usage()
{
  echo "Usage: copyright.sh [-h|-help] [-options]"
}

FullUsage()
{
  Usage
  echo "  options:"
  echo "  -h|-help        this help screen"
  echo "  -check_exe      checks executable files for "
  echo "  -check_tex      this help screen"
  echo "  -check_special  this help screen"
  echo "... options not finished..."
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

show_copy="NO"
find_type="none"
write="NO"

while [ -n "$1" ]
do
   case $1 in
        -h|-help)       FullUsage; exit 0;;
        -check_exe)     what="check_exe";;
        -check_tex)     what="check_tex";;
        -check_special) what="check_special";;
        -check_str)     what="check_str";;
        -check_fortran) what="check_fortran";;
        -check_c)       what="check_c";;
        -check_all)     what="check_all";;
        -find_type)     what="find_type"; find_type=$2; shift;;
        -print_type)    what="print_type";;
        -show_copy)     what="show_copy";;
        -show_stats)    what="show_stats";;
        -write)         write="YES";;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   shift
done

if [ -n "$1" ]; then	#extra argument
  echo "no extra argument allowed: $1"
  Usage; exit 0
elif [ -z "$what" ]; then
  Usage; exit 0
fi

echo "running in directory: $PWD"

#---------------------------------------------------------------

if [ -z "$what" ]; then
  Usage; exit 0
elif [ $what = "print_type" ]; then
  PrintFileType
elif [ $what = "find_type" ]; then
  FindFileType
elif [ $what = "check_exe" ]; then
  CheckExeType
elif [ $what = "check_tex" ]; then
  CheckTexType
elif [ $what = "check_str" ]; then
  CheckStrType
elif [ $what = "check_fortran" ]; then
  CheckFortranType
elif [ $what = "check_c" ]; then
  CheckCType
elif [ $what = "check_special" ]; then
  CheckSpecialType
elif [ $what = "show_copy" ]; then
  ShowCopyright
elif [ $what = "check_all" ]; then
  CheckAllType
elif [ $what = "show_stats" ]; then
  ShowStats
else
  Usage
fi

#---------------------------------------------------------------

