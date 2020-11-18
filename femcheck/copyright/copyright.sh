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
# external routines used:
#
#	revision_log.sh			all programs on handling revision log
#	revision_log.pl			all programs on handling revision log
#	extension_stats.pl		compute statistics on extensions
#	find_file_type.pl		find file type of file
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

  files=$( findf -x '*' | grep -v '/tmp/' | grep -v '/arc/' \
		| grep -v /.git/ )
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

  files=$( findf '*' | grep -v '/tmp/' | grep -v '/arc/' \
		| grep -v /.git/ )

  for file in $files
  do
    [ -d $file ] && continue
    IsBinary $file
    [ $? -ne 0 ] && echo "skipping binary file: $file " && continue
    #if IsBinary $file; then echo "skipping binary file: $file "; continue; fi
    first=$( head -1 $file )
    if [[ $first =~ '#!/'.* ]]; then
      if [ ! -x $file ]; then
	errors=$(( errors + 1 ))
        if [ $write = "YES" ]; then
          echo "making $file executable"
	  chmod +x $file
        else
          echo "*** $file is a script but is not executable"
        fi
      fi
    fi
  done
  if [ $errors -ne 0 ]; then
    if [ $write = "YES" ]; then
      echo "$errors files have been made executable..."
      errors=0
    else
      echo "use -write to make scripts executable"
      exit 1
    fi
  fi

  echo "--- printing scripts that have no copyright"

  files=$( findf '*' | grep -v '/tmp/' | grep -v '/arc/' \
		| grep -v /.git/ | grep -v /GD/ \
		| grep -v /Mail-Sender | grep -v /femersem/ )

  HandleCopyright script
}

CheckTexType()
{
  echo "================================================"
  echo "--- CheckTexType: checking tex files"
  echo "================================================"

  echo "--- printing all file types that are tex files"

  files=""
  MakeFilesFromExt tex
  FilterFiles /tmp/ /arc/ 

  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright tex
}

CheckStrType()
{
  echo "================================================"
  echo "--- CheckStrType: checking str files"
  echo "================================================"

  echo "--- printing all file types that are str files"

  files=""
  MakeFilesFromExt str
  FilterFiles /tmp/ /arc/ 

  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright text
}

CheckCType()
{
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

  HandleCopyright c
}

CheckFortranType()
{
  echo "================================================"
  echo "--- CheckFortranType: checking fortran files"
  echo "================================================"

  echo "--- printing all file types that are fortran files"

  files=""
  MakeFilesFromExt f F90 f90 F h
  FilterFiles /tmp/ /arc/ 
  #FilterDirs femgotm '.*.F90' '.*.h'
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

  HandleCopyright fortran
}

CheckSpecialType()
{
  echo "================================================"
  echo "--- CheckSpecialType: checking special files"
  echo "================================================"

  echo "--- printing all file types that are special files"

  special="Makefile makefile README Rules.make Include.make"

  MakeFilesFromName $special
  FilterFiles /tmp/ /arc/ 
  FilterFiles /Mail-Sender /femersem/

  if [ "$files" != "" ]; then
    echo $files \
	| xargs file -b \
	| sort -u
  fi

  echo "--- printing files that have no copyright"

  HandleCopyright text
}

CheckAllType()
{
  CheckExeType
  CheckTexType
  CheckStrType
  CheckCType
  CheckFortranType
  CheckSpecialType

  files=""
  CheckRev $extra
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

MakeFilesFromCommandLine()
{
  if [ -z "$files" ]; then
    files=$( findf '*' )
  elif [[ $files == .* ]]; then	#starts with ., therefore is extension
    exts=$files
    files=""
    MakeFilesFromExt $exts
  fi
}

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
    ext=$( echo $ext | sed -e 's/^\.//' )	#eliminate dot
    aux=$( findf '*.'$ext )
    files="$files $aux"
  done
}

FilterExtensions()
{
  local filtered=""

  for ext
  do
    ext=$( echo $ext | sed -e 's/^\.//' )	#eliminate dot
    ext="\\.$ext"
    aux=$( echo $files | tr ' ' '\n' | grep "$ext$" )
    filtered="$filtered $aux"
    aux=$( echo $files | tr ' ' '\n' | grep -v "$ext$" )
    files=$aux
  done

  #filtered=$( echo $filtered | tr '\n' ' ' )
  #echo "filtered files: $filtered"
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

IsBinary()
{
  local file=$1

  local filename=$(basename -- "$file")
  local extension="${filename##*.}"

  [ $extension = "pdf" ] && return 0
  [ $filename = "CR" ] && return 0

  grep -qI . $1

  local status=$?
  if [ $status -ne 0 ]; then
    :
    #echo "  *** file might be binary: $file"
  else
    :
  fi
  return $status
}

GetFileType()
{
  $copydir/find_file_type.pl $1
}

DetermineFileType()
{
  MakeFilesFromCommandLine
  FilterFiles /tmp/ /arc/ /orig/
  FilterFiles /.git/
  FilterFiles /femersem/

  for file in $files
  do
    [ -d $file ] && continue
    [ -L $file ] && continue
    type=$( $copydir/find_file_type.pl $file )
    echo "$file: $type"
  done
}

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
  type=$1
  errors=0

  findtype=NO
  [ -z "$type" ] && findtype=YES

  for file in $files
  do
    [ -d $file ] && continue
    IsBinary $file
    [ $? -ne 0 ] && echo "skipping binary file: $file " && continue
    [ $findtype = YES ] && type=$( GetFileType $file )
    if [ "$type" = "script" ]; then
      first=$( head -1 $file )
      [[ ! $first =~ '#!/'.* ]] && continue
    fi

    $copydir/revision_log.pl -onlycopy $file
    error=$?
    if [ $error -ne 0 ]; then
      if [ $write = "YES" ]; then
        echo "*** $file has no copyright... inserting"
        $copydir/revision_log.sh -newcopy -write $file
      else
        echo "*** $file has no copyright..."
        errors=$(( errors + 1 ))
      fi
    fi
  done

  if [ $errors -gt 0 -a $write = "NO" ]; then
    echo "$errors files have error in copyright... use -write to insert"
    exit 1
  fi
}

DoCopyright()
{
  if [ $# -eq 0 ]; then
    action=show
  elif [ $1 = -show ]; then
    action=show
  elif [ $1 = -check ]; then
    action=check
  else
    echo "ShowCopyright: no such action: $1"
    return
  fi

  errors=0
  show_copy="YES"

  #echo "files: $files"
  MakeFilesFromCommandLine
  FilterFiles /tmp/ /arc/ /orig/
  FilterFiles /.git/ /femersem/
  FilterFiles /INPUT/
  FilterFiles /oceanlib/ /oceanlib_txt/
  FilterFiles /Mail-Sender-0.8.13/ /codepage/ /GD
  FilterExtensions pdf ps eps jpg gif png
  FilterExtensions dat grd

  for file in $files
  do
    [ -d $file ] && continue
      error=0
      head -50 $file | grep -E "^..\s*This file is part of SHYFEM." > /dev/null
      [ $? -ne 0 ] && error=$(( error + 1 ))
      head -50 $file | grep -E "^..\s*Copyright \(C\)" > /dev/null
      [ $? -ne 0 ] && error=$(( error + 10 ))
      #echo "-------------- $file $error"
      if [ $error -eq 0 ]; then
        if [ $action = show ]; then
          echo "$file has copyright"
        fi
      else
        if [ $action = check ]; then
          echo "$file has no copyright"
	  newfiles="$newfiles $file"
        fi
      fi
  done

  [ -z "$newfiles" ] && return
  [ $write = NO ] && return

  files=$newfiles
  echo "integrating copyright in files"
  HandleCopyright
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

  $copydir/extension_stats.pl $files
}

CheckRev()
{
  option=$extra
  [ -z "$option" ] && option="-check"
  
  echo "================================================"
  echo "--- CheckRev: checking revision log"
  echo "================================================"

  MakeFilesFromCommandLine
  FilterFiles /tmp/ /arc/ /orig/
  FilterFiles /.git/
  FilterFiles /oceanlib/ /oceanlib_txt/
  FilterFiles /Mail-Sender-0.8.13/ /codepage/ /GD
  FilterFiles /femersem/
  FilterFiles /bugs/ /INPUT/
  FilterFiles /copyright_ .swp

  $copydir/revision_log.sh $option $files
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
  echo "Usage: copyright.sh [-h|-help] [-options] [files]"
}

FullUsage()
{
  Usage
  echo "  options:"
  echo "  -h|-help           this help screen"
  echo "  -check_exe         checks executable files for coherence"
  echo "  -check_tex         checks tex files for coherence"
  echo "  -check_special     checks special files for coherence"
  echo "  -check_str         checks str files for coherence"
  echo "  -check_fortran     checks fortran files for coherence"
  echo "  -check_c           checks c files for coherence"
  echo "  -check_all         checks all files for coherence"
  echo "  -find_type type    finds files with type type"
  echo "  -print_type        prints types of files"
  echo "  -determine_type    determines types of files"
  echo "  -check_copy        checks if files have copyright notice"
  echo "  -show_copy         shows if files have copyright notice"
  echo "  -show_stats        shows statistics of file extensions"
  echo "  -check_rev         checks revision log"
  echo "  -write             if missing, insert copyright"
  echo "  files can be list of files (*.f) or extension (.f)"
  echo "  with no files all files are searched, search is recursive"
  echo "  extra options for -check_rev;"
  echo "  --check            check revision log (default)"
  echo "  --gitrev           integrates git revision log into file"
  echo "  --gitmerge         merges git revision log into file"
  echo "  --stats            show stats of revision log and copyright"
  #echo "  --crewrite         rewrites c style revision log to standard"
  echo "  --substdev         substitutes developer name with new name"
  echo "  --onlycopy         only checks copyright notice"
  echo "  --updatecopy       updates copyright with info from revision log"
  echo "  --newcopy          substitute old with new copyright"
  echo "  --gui              if files are changed show in gui diff"
  echo "  --write            really write changes to file"
  echo "  --keep             keep changed files (.revnew) for inspection"
  echo "  for general check use -check_rev, -check_all, -check_copy"
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
        -h|-help)         FullUsage; exit 0;;
        -check_exe)       what="check_exe";;
        -check_tex)       what="check_tex";;
        -check_special)   what="check_special";;
        -check_str)       what="check_str";;
        -check_fortran)   what="check_fortran";;
        -check_c)         what="check_c";;
        -check_all)       what="check_all";;
        -find_type)       what="find_type"; find_type=$2; shift;;
        -print_type)      what="print_type";;
        -determine_type)  what="determine_type";;
        -check_copy)      what="check_copy";;
        -show_copy)       what="show_copy";;
        -show_stats)      what="show_stats";;
        -check_rev)       what="check_rev";;
        -write)           write="YES";;
        --*)              extra="$extra ${1#?}";;		#pop one -
        -*)               ErrorOption $1; exit 1;;
        *)                break;;
   esac
   shift
done

if [ -n "$1" ]; then	#extra argument
  files=$*
elif [ -z "$what" ]; then
  Usage; exit 0
fi

echo "looking for files in directory: $PWD"

#---------------------------------------------------------------

if [ -z "$what" ]; then
  Usage; exit 0
elif [ $what = "print_type" ]; then
  PrintFileType
elif [ $what = "find_type" ]; then
  FindFileType
elif [ $what = "determine_type" ]; then
  DetermineFileType
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
  DoCopyright -show
elif [ $what = "check_copy" ]; then
  DoCopyright -check
elif [ $what = "check_all" ]; then
  CheckAllType
elif [ $what = "show_stats" ]; then
  ShowStats
elif [ $what = "check_rev" ]; then
  CheckRev $extra
else
  Usage
fi

#---------------------------------------------------------------

