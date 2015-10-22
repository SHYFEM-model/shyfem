#!/bin/sh

log=CHECKLOG
rm -f $log

missing=""

#---------------------------------------------------

CheckCommand()
{
  name=$1
  command=$2
  error=$3
  options=$4

  if [ -z "$error" ]; then
    error=0
  fi

  ($command) >> $log 2>&1

  status=$?

  if [ $status -eq $error ]; then
    echo "... $name is ${green}installed${normal}"
  else
    if [ "$options" != "quiet" ]; then
      echo "*** $name is ${red}not installed${normal}"
    fi
    missing="$missing $name"
  fi
}

CreateInputFiles()
{

cat > test.f <<EOI
	include 'netcdf.inc'
        write(6,*) 'Hello, world.'
        end
EOI

cat > test.c <<EOI
#include <stdio.h>
#include <X11/Xlib.h>

int main( void )
{
  printf("Hello world.\n");
  return 0;
}
EOI

echo "quit" > quit.tmp
}

CleanUp()
{
	[ -f test.f ] && rm -f test.f
	[ -f test.c ] && rm -f test.c
	[ -f quit.tmp ] && rm -f quit.tmp
	[ -f a.out ] && rm -f a.out
}

#---------------------------------------------------

CheckFortranCompiler()
{
  #local fortran_available=""
  local missing_save=$missing

  CheckCommand g77 "g77 -v" "" "quiet"
  [ $status -eq 0 ] && fortran_available="$fortran_available g77"
  CheckCommand gfortran "gfortran -v" "" "quiet"
  [ $status -eq 0 ] && fortran_available="$fortran_available gfortran"
  CheckCommand ifort "ifort -v" "" "quiet"
  [ $status -eq 0 ] && fortran_available="$fortran_available ifort"

  #fortran_available=""		# fake error

  if [ -n "$fortran_available" ]; then
    echo "... the following Fortran compilers are available:"
    echo "          $fortran_available"
    echo "... please set COMPILER in Rules.make to one of the compilers above"
  else
    echo "*** No Fortran compiler found... please install a fortran compiler"
    echo "    (on Linux you can easily install gfortran)"
    missing_save="$missing_save g77 gfortran ifort"
  fi

  fortran=`echo $fortran_available | sed -e 's/ .*//'`	#first fortran found
  missing=$missing_save
}

CheckX11()
{
  CheckCommand "X11-develop" \
		"gcc -L/usr/X11/lib -L/usr/X11R6/lib -lXt -lX11 test.c"

  if [ $status -ne 0 ]; then
    echo "*** X11 development package must be installed"
    echo "    please try with: libx11 libx11-common libx11-dev libxt-dev"
  fi
}

CheckNetcdf()
{
  netcdf=`GetMacro NETCDF`
  netcdfdir=`GetMacro NETCDFDIR`

  #echo "netcdf: $netcdf $netcdfdir    $fortran"

  [ -z "$fortran" ] && return		# no fortran compiler
  [ "$netcdf" != "true" ] && return	# no netcdf requested

  CheckCommand netcdf "$fortran -L$netcdfdir/lib -I$netcdfdir/include test.f"

  if [ $status -ne 0 ]; then
    echo "*** netcdf seems not to be installed"
    echo "    If you need netcdf please install libnetcdf-dev or similar"
  fi
}

GetMacro()	# gets macro from Rules.make file
{
  what=$1
  rules="../Rules.make"
  [ $# -gt 1 ] && rules="$2"

  macro=$( cat $rules | sed -e 's/ *//g' | grep "^$what=" \
	| tail -1 | sed -e 's/.*=//')

  echo "$macro"
}

SetColors()
{
  ncolors=$(tput colors)
  if test -n "$ncolors" && test $ncolors -ge 8; then
        bold="$(tput bold)"
        underline="$(tput smul)"
        standout="$(tput smso)"
        normal="$(tput sgr0)"
        black="$(tput setaf 0)"
        red="$(tput setaf 1)"
        green="$(tput setaf 2)"
        yellow="$(tput setaf 3)"
        blue="$(tput setaf 4)"
        magenta="$(tput setaf 5)"
        cyan="$(tput setaf 6)"
        white="$(tput setaf 7)"
  fi
}

#---------------------------------------------------

SetColors
CreateInputFiles

echo
echo "... ${bold}checking base applications... (needed)${normal}"

CheckCommand make "make -v"
CheckCommand bash "bash --version"
CheckCommand perl "perl -v"

echo
echo "... ${bold}checking Fortran compilers (needed)${normal}"

CheckFortranCompiler

echo
echo "... ${bold}checking c compiler and X11 (needed)${normal}"

CheckCommand gcc "gcc -v"
CheckCommand g++ "g++ -v"
CheckX11

echo
echo "... ${bold}commodity programs (recommended)${normal}"
CheckCommand dialog  "dialog --version"

echo
echo "... ${bold}checking graphical routines (recommended)${normal}"

CheckCommand "ghostview (gv)" "gv --version"
CheckCommand "ghostscript (gs)" "gs -v"
CheckCommand gnuplot "gnuplot quit.tmp"
CheckCommand "ImageMagick (imagemagick)" "mogrify -version"
CheckCommand gifsicle "gifsicle --version"
CheckNetcdf

echo
echo "... ${bold}checking additional routines (not urgently needed)${normal}"

CheckCommand "latex (texlive-full)" "latex -v"
CheckCommand dvips "dvips -v"
CheckCommand python "python -V"

CheckCommand ssh "ssh -V"
CheckCommand at "at -l"

#CheckCommand dnotify "dnotify --version"
#CheckCommand "Acrobat Reader" "acroread -help"

if [ -n "$missing" ]; then
  echo ""
  echo "${bold}The following programs seem not be installed: ${normal}"
  echo ""
  echo    "${red}$missing${normal}"
  echo ""
  echo "Please see the messages above to find out if these programs"
  echo "are indispensable to run the model."
  echo "It might be a good idea to install these programs anyway."
  echo "How exactly this is done depends a lot on the type of"
  echo "distribution you are using. You could use a software management"
  echo "tool that comes with your distribution."
  echo "Search for the keyword of the missing file and then install"
  echo "the files on your harddisk."
  echo ""
else
  echo ""
  echo    "${green}All programs installed${normal}"
  echo ""
fi

CleanUp

# xanim

