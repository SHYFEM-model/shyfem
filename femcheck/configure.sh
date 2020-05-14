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
# configures Rules.make file
#
#-----------------------------------------------------

TEMP=tmp.tmp

femdir=..
rules_make_file="$femdir/Rules.make"
configurepl="$femdir/femcheck/configure.pl"

fcompilers="GNU_G77 GNU_GFORTRAN INTEL IBM PORTLAND"
ccompilers="GNU_GCC INTEL IBM"
solvers="GAUSS SPARSKIT"
ecos="NONE EUTRO ERSEM AQUABC"

#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

create_fortran_list()
{
  list=""
  for comp in $fcompilers
  do
    case $comp in
    	GNU_G77 )	include_in_list $comp "g77 --version" ;;
    	GNU_GFORTRAN )	include_in_list $comp "gfortran --version" ;;
    	INTEL )		include_in_list $comp "ifort --version" ;;
    	IBM )		include_in_list $comp "xlf_r" ;;
    	PORTLAND )	include_in_list $comp "pgf90 -V" ;;
	* )		raise_error "$comp not in list" ;;
    esac
  done
  #show_message_time 2 "list: $list"
  fcompilers=$list

  has_substring $fcompiler $fcompilers
  [ $result = false ] && fcompiler=NONE
}

create_c_list()
{
  list=""
  for comp in $ccompilers
  do
    case $comp in
    	GNU_GCC )	include_in_list $comp "gcc --version" ;;
    	INTEL )		include_in_list $comp "icc --version" ;;
    	IBM )		include_in_list $comp "xlc" ;;
	* )		raise_error "$comp not in list" ;;
    esac
  done
  #show_message_time 2 "list: $list"
  ccompilers=$list

  has_substring $ccompiler $ccompilers
  [ $result = false ] && ccompiler=NONE
}

read_rules_make()
{
  if [ ! -f $rules_make_file ]; then
    echo "Cannot find file: $rules_make_file"
    exit 1
  fi

  fcompiler=$( read_macro FORTRAN_COMPILER )
  ccompiler=$( read_macro C_COMPILER )
  parallel_omp=$( read_macro PARALLEL_OMP )
  solver=$( read_macro SOLVER )
  netcdf=$( read_macro NETCDF )
  netcdfdir=$( read_macro NETCDFDIR )
  gotm=$( read_macro GOTM )
  eco=$( read_macro ECOLOGICAL )
}

read_macro()
{
  result="UNKNOWN"
  line=$( grep "^\s*$1\s*=" $rules_make_file | tail -1 )
  macro=$( echo $line | sed -E 's/.*=\s*//' )
  macro=$( echo $macro | sed -E 's/\s*\#.//' )
  #echo "$line ... $macro"
  echo "$macro"
}

write_rules_make()
{
  
  if [ $fcompiler = NONE ]; then
    raise_error "You must choose a fortran compiler"
  fi
  if [ $ccompiler = NONE ]; then
    raise_error "You must choose a c compiler"
  fi

  show_message "writing Rules.make file"
  #show_message "writing Rules.make file\n...not yet ready"

  $configurepl \
		-FORTRAN_COMPILER=$fcompiler \
		-C_COMPILER=$ccompiler \
		-PARALLEL_OMP=$parallel_omp \
		-SOLVER=$solver \
		-NETCDF=$netcdf \
		-NETCDFDIR=$netcdfdir \
		-GOTM=$gotm \
		-ECOLOGICAL=$eco \
		$rules_make_file > $TEMP

  cp --backup=numbered $TEMP $rules_make_file
}

check_needed_software()
{
  line=""
  error=NO

  if [ -z "$fcompilers" ]; then
    show_message_ok "Cannot find a fortran compiler"
    error=YES
  fi

  if [ -z "$ccompilers" ]; then
    show_message_ok "Cannot find a c compiler"
    error=YES
  fi

  if [ error = YES ]; then
    show_message_time 2 "Please install and re-run configure"
    clean_up
  fi

  is_command_installed "bash --version"
  [ $result = false ] && line="$line\nbash"

  is_command_installed "make --version"
  [ $result = false ] && line="$line\nmake"

  is_command_installed "perl --version"
  [ $result = false ] && line="$line\nperl"

  is_command_installed "g++ --version"
  [ $result = false ] && line="$line\nperl"

  if [ -n "$line" ]; then
    show_message_ok "Some necessary programs are not installed:\n$line"
    error=YES
  fi

  if [ error = YES ]; then
    show_message_time 2 "Please install and re-run configure"
    clean_up
  fi

  is_lib_installed libX11
  [ $result = false ] && line="$line\nlibX11"

  is_lib_installed libXt
  [ $result = false ] && line="$line\nlibXt"

  if [ -n "$line" ]; then
    show_message_ok "Some recommended libraries are not installed:\n$line"
  fi
}

check_dialog()
{
  is_command_installed "dialog --version"
  if [ $result = false ]; then
    echo 'to run "make configure" the program "dialog" is needed'
    echo 'please install "dialog" and then re-run "make configure"'
    #sleep 2
    exit 0
  fi
}

#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

process_netcdf()
{
  is_lib_installed libnetcdff
  if [ $result = false ]; then
    $netcdf=false
    show_message_time 2 "netcdf library seems not to be installed"
    return
  fi

  get_lib_dir libnetcdff
  netcdfdir=$result

  process_yesno "Use netcdf library in $netcdfdir?"
  netcdf=$result
}

process_gotm()
{
  process_yesno "Use GOTM library?"
  gotm=$result
}

process_parallel()
{
  process_yesno "Compile for parallel execution?"
  parallel_omp=$result
}

process_fortran()
{
  process_list "Choose Fortran Compiler:" $fcompiler $fcompilers
  fcompiler=$result
}

process_c()
{
  process_list "Choose C Compiler:" $ccompiler $ccompilers
  ccompiler=$result
}

process_eco()
{
  process_list "Choose Ecological Model:" $eco $ecos
  eco=$result
}

process_solver()
{
  process_list "Choose Matrix Solver:" $solver $solvers
  solver=$result
}

not_ready()
{
  show_message "Not ready yet..."
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

process_list()
{
  message=$1
  defvalue=$2
  shift; shift
  num=$#
  line=""
  array=("${@}")

  for i in $(seq 1 $num)
  do
    index=$((i-1))
    entry=${array[$index]}
    status=$(on_off $entry $defvalue)
    echo "$i $index $status"
    line="$line $i $entry $status"
  done

  #echo $line; exit 1

  dialog --title "$message" \
    --radiolist "$message" 11 60 $num \
    $line \
    2>$TEMP

  i=`cat $TEMP`
  index=$((i-1))
  #show_message "i returned $i   index $index"

  result=$defvalue
  [ $index -eq -1 ] && return
  result=${array[$index]}
}

process_yesno()
{
  dialog --title "Message"  --yesno "$1" 6 25
  if [ $? -eq 0 ]; then
    result="true"
  else
    result="false"
  fi
  #echo $result
}

show_message_ok()
{
  dialog --msgbox "$1" 10 30
}

show_message()
{
  show_message_time 1 "$1"
}

show_message_time()
{
  dialog --infobox "$2" 10 30 ; sleep $1
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

clean_up() {		 # clean up and exit
  clear
  rm -f $TEMP
  exit
}

on_off() { 		# utility function
  if [ "$1" = "$2" ] ; then echo on ; else echo off ; fi
}

get_lib_dir()
{
  line=$( ldconfig -p | grep $1 | tail -1 )
  dir=$( echo $line | sed -E 's/^.*\s//' | sed -E 's/\/lib.*//' )
  result=$dir
}

is_lib_installed()
{
  ldaux=ldconfig
  $ldaux > /dev/null 2>&1
  [ $? -eq 127 ] && ldaux=/sbin/ldconfig

  n=$( $ldaux -p | grep $1 | wc -l )
  result=false
  [ $n -ne 0 ] && result=true
}

is_command_installed()
{
  $1 > /dev/null 2>&1
  if [ $? -eq 127 ]; then
    result=false
  else
    result=true
  fi
}

has_substring()
{
  string=$1
  list=$2

  if echo "$list" | grep -q "$string"; then
    result=true
  else
    result=false
  fi
}

include_in_list()
{
  entry=$1
  command=$2
  is_command_installed "$command"
  if [ $result = true ]; then
    list="$list $entry"
  fi
}

raise_error()
{
  show_message_ok "ERROR: $1"
  clean_up
}

#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

check_dialog
read_rules_make
create_fortran_list
create_c_list
check_needed_software

while [ 1 ]
do
  dialog --menu "Choose entry to change:" 18 50 9 \
		1 "Fortran compiler:    $fcompiler" \
		2 "C compiler:          $ccompiler" \
		3 "Parallel:            $parallel_omp" \
		4 "Matrix solver:       $solver" \
		5 "Use netcdf:          $netcdf" \
		6 "GOTM library:        $gotm" \
		7 "Ecological model:    $eco" \
		8 "quit without saving" \
		9 "save and exit" \
		2>$TEMP
  choice=`cat $TEMP`
  #show_message "Your choice was: $choice"
  case $choice in
	    	1 )		process_fortran ;;
	    	2 )		process_c ;;
	    	3 )		process_parallel ;;
	    	4 )		process_solver ;;
    		5 )		process_netcdf ;;
    		6 )		process_gotm ;;
	    	7 )		process_eco ;;
	    	8 )		break ;;
	    	9 )		write_rules_make; break ;;
	    	* )		break ;;
  esac
done

clean_up

#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

