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
# command line options - use as module in other bash scripts
#
# Copyright (c) 2014-2015, 2020 Georg Umgiesser - ISMAR-CNR
#
# version:
#
# 1.00	26.01.2015	routine written from scratch
# 1.10	24.03.2015	new routine clo_shift_item() and clo_get_shifted()
# 1.20	17.11.2020	handle spaces in item
#
# please see clo_test at end of file for usage information
# please insert ". clo.sh" close to beginning of the bash script
#
#------------------------------------------------------

# global variables

declare -A clo_option		# value for option
declare -A clo_description	# description for option
declare -A clo_expect		# expect extra argument
declare -A clo_exists		# if YES then option exists
declare -a clo_info		# information before single options
declare -a clo_extra		# extras info after single options
declare -a clo_sequence		# sequence in which options have been added

declare clo_items=""		# items after options
declare clo_shifted=""		# after shift contains shifted argument
declare clo_lenmax=11		# max lenght of name and expect (-v|-version)
declare clo_routine=""		# routine name
declare clo_item_info=""	# item info for usage
declare clo_version=""		# version of routine

declare clo_debug="NO"		# debug?

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------

clo_init()
{
  clo_routine=$(basename $0)
  clo_item_info=$1
  clo_version=$2
}

clo_parse_options()
{
  local name exists expect
  [ $clo_debug = "YES" ] && echo "parsing options..."

  clo_items=("${@}")

  while [ -n "${clo_items[0]}" ]
  do
     item="${clo_items[0]}"
     name=`echo $item | sed -E 's|^-||'`		# strip - from option
     [ "$name" = "$item" ] && break			# no option
     exists=${clo_exists[$name]}
     if [ "$name" = "help" -o "$name" = "h" ]; then
       clo_fullusage
       exit 1
     elif [ "$name" = "version" -o "$name" = "v" ]; then
       clo_usage
       echo "version = $clo_version"
       exit 1
     elif [ "$exists" != "YES" ]; then
       echo "Unknown option: $1"
       clo_usage
       exit 1
     fi
     expect=${clo_expect[$name]}
     if [ -n "$expect" ]; then
       clo_option[$name]=$2
       shift
     else
       clo_option[$name]="YES"
     fi
     [ $clo_debug = "YES" ] && echo "option found: $name = ${clo_option[$name]}"
     shift
     clo_items=("${clo_items[@]:1}")
     clo_export_option $name
  done

  [ $clo_debug = "YES" ] && echo "after parsing options: $*"
}

clo_usage()
{
  echo "Usage: $clo_routine [-h|-help] [-options] $clo_item_info"
}

clo_fullusage()
{
  local name i print expect descr
  clo_usage
  for ((i=0;i<${#clo_info[*]};i++)); do
    echo "  ${clo_info[$i]}"
  done
  echo "  options:"
  clo_print_option "-h|-help" "this help screen"
  clo_print_option "-v|-version" "print version of routine"
  for ((i=0;i<${#clo_sequence[*]};i++)); do
    name=${clo_sequence[$i]}
    print=$name
    expect=${clo_expect[$name]}
    [ -n "$expect" ] && print="$name $expect"
    descr=${clo_description[$name]}
    clo_print_option "-$print" "$descr"
  done
  for ((i=0;i<${#clo_extra[*]};i++)); do
    echo "  ${clo_extra[$i]}"
  done
}

clo_print_option()
{
  local text descr m print
  text=$1
  descr=$2
  m=$((clo_lenmax+3))
  print=$(printf "%-${m}s" "$text")
  echo "  $print $descr"
}

#------------------------------------------------------

clo_add_option()
{
  local aux name expect a val descr
  if [ $# -ne 3 ]; then
    echo "*** Usage: clo_add_option var value description"
    exit 1
  fi
  aux=$1
  [ ${#aux} -gt $clo_lenmax ] && clo_lenmax=${#aux}
  name=""
  expect=""
  for a in $aux
  do
    if [ "$name" = "" ]; then
      name=$a
    else
      expect=$a
    fi
  done
  val=$2
  descr=$3
  clo_option[$name]=$val
  clo_description[$name]=$descr
  clo_expect[$name]=$expect
  clo_exists[$name]="YES"
  clo_sequence=("${clo_sequence[@]}" "$name")
  clo_export_option $name
}

clo_export_option()
{
  local opt val
  if [ $# -eq 0 ]; then
    echo "*** Usage: clo_set_option var"
    exit 1
  fi
  opt=$1
  val=${clo_option[$opt]}
  eval $opt=$val
}

clo_get_option()
{
  local opt val
  if [ $# -eq 0 ]; then
    echo "*** Usage: clo_get_option var"
    exit 1
  fi
  opt=$1
  val=${clo_option[$opt]}
  echo $val
}

clo_show_info()
{
  local text name i val expect descr
  text="show info $1"
  echo "--------------------------------------"
  echo "$text"
  echo "--------------------------------------"
  echo "available options:"
  for ((i=0;i<${#clo_sequence[*]};i++)); do
    name=${clo_sequence[$i]}
    val=${clo_option[$name]}
    expect=${clo_expect[$name]}
    [ -n "$expect" ] && expect="($expect)"
    descr=${clo_description[$name]}
    echo "  $name = $val  descr = $descr   $expect"
  done
  echo "available infos:"
  for ((i=0;i<${#clo_info[*]};i++)); do
    echo "  $i : ${clo_info[$i]}"
  done
  echo "available extras:"
  for ((i=0;i<${#clo_extra[*]};i++)); do
    echo "  $i : ${clo_extra[$i]}"
  done
  echo "--------------------------------------"
}

clo_add_extra()
{
  clo_extra=("${clo_extra[@]}" "$1")
}

clo_add_info()
{
  clo_info=("${clo_extra[@]}" "$1")
}

#------------------------------------------------------

clo_number_of_items()
{
  echo ${#clo_items[@]}
}

clo_print_items()
{
  echo "${clo_items[@]}"
}

clo_shift_item()
{
  local shifted=""
  local n_items=${#clo_items[@]}

  if [ $n_items -gt 0 ]; then
    clo_shifted=${clo_items[0]}
    clo_items=("${clo_items[@]:1}")
  fi
}

clo_get_shifted()
{
  echo "$clo_shifted"
}

clo_check_items()
{
  local n_items=${#clo_items[@]}

  if [ $# -ne 1 ]; then
    echo "*** Usage: clo_check_items min-number-of-items"
    exit 1
  elif [ $n_items -lt "$1" ]; then
    #echo "Expecting at least $1 items"
    clo_usage
    exit 1
  fi
}

#------------------------------------------------------

clo_print_debug()
{
  echo "options:"
  for key in "${!clo_exists[@]}"
  do 
    val=${clo_option[$key]}
    descr=${clo_description[$key]}
    expect=${clo_expect[$key]}
    echo "$key - $val - $descr - $expect"
  done

  echo "items:"
  for i in ${!clo_items[@]}
  do
    echo "$i - ${clo_items[$i]}"
  done
}

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------

# fortran interface:
#        call clo_add_info('returns info on a fem file')
#        call clo_add_option('write',.false.,'write min/max of values')
#        call clo_add_option('quiet',.false.,'do not be verbose')
#        call clo_add_option('itmin time',"-1",'do not be verbose')
#        call clo_parse_options(1)  !expecting (at least) 1 file after options
#        call clo_get_option('write',bwrite)
#        call clo_get_option('quiet',bquiet)
#        nfile = clo_number_of_files()

clo_test()
{
  clo_init "fem-file(s)" "1.0"
  clo_add_info "returns info on a fem file"
  clo_add_option write "NO" "write min/max of values"
  clo_add_option quiet "NO" "do not be verbose"
  clo_add_option debug "NO" "debug mode"
  clo_add_option "itmin time" "-1" "only process starting from time"
  clo_add_option "itmax time" "-1" "only process until time"
  clo_add_extra "format for time is YYYYMMDD or 'YYYY-MM-DD::hh:mm:ss'"
  clo_add_extra "also relative time can be given (in seconds)"
  #clo_show_info "before parsing"
  clo_parse_options "$@"
  clo_check_items 1
  [ "$debug" = "YES" ] && clo_show_info "after parsing"
  items=$(clo_print_items)
  for item in $items
  do
    echo "processing item: $item"
  done
}

clo_main()
{
  calling_prog=`echo $0 | sed -E "s/.*\///"`
  [ "$calling_prog" = "clo.sh" ] && clo_test "$@"
}

#------------------------------------------------------

clo_main "$@"

#------------------------------------------------------
