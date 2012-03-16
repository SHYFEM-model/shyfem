#!/bin/sh

log=CHECKLOG
rm -f $log

missing=""

CheckCommand()
{
  name=$1
  command=$2
  error=$3

  if [ -z "$error" ]; then
    error=0
  else
    echo "error passed ............"
  fi

  ($command) >> $log 2>&1

  status=$?

  if [ $status -eq $error ]; then
    echo "... $name is installed"
  else
    echo "*** $name is not installed"
    missing="$missing $name"
  fi
}

CreateInputFiles()
{

cat > test.f <<EOI
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

#---------------------------------------------------

CheckFortranCompiler()
{
  local fortran_available=""
  local missing_save=$missing

  CheckCommand g77 "g77 -v"
  [ $status -eq 0 ] && fortran_available="$fortran_available g77"
  CheckCommand gfortran "gfortran -v"
  [ $status -eq 0 ] && fortran_available="$fortran_available gfortran"
  CheckCommand ifort "ifort -v"
  [ $status -eq 0 ] && fortran_available="$fortran_available ifort"

  #fortran_available=""		# fake error

  if [ -n "$fortran_available" ]; then
    echo "... the following Fortran compilers are available:"
    echo "          $fortran_available"
    echo "... please set COMPILER in Rules.make to one of the compilers above"
  else
    echo "*** No Fortran compiler found... please install a fortran compiler"
    missing_save="$missing_save g77 gfortran ifort"
  fi

  missing=$missing_save
}

CheckX11()
{
  CheckCommand X11 "gcc -L/usr/X11/lib -L/usr/X11R6/lib -lXt -lX11 test.c"

  if [ $status -ne 0 ]; then
    echo "*** X11 development package must be installed"
    echo "    please try with: libx11 libx11-common libx11-dev libxt-dev"
  fi
}

#---------------------------------------------------

CreateInputFiles

echo
echo "... checking base applications... (needed)"

CheckCommand make "make -v"
CheckCommand bash "bash --version"
CheckCommand perl "perl -v"

echo
echo "... checking Fortran compilers (needed)"

CheckFortranCompiler

echo
echo "... checking c compiler and X11 (needed)"

CheckCommand gcc "gcc -v"
CheckCommand g++ "g++ -v"
CheckX11

echo
echo "... checking graphical routines (recommended)"

CheckCommand "ghostview (gv)" "gv --version"
CheckCommand ghostscript "gs -v"
CheckCommand gnuplot "gnuplot quit.tmp"
CheckCommand ImageMagic "mogrify -version"
CheckCommand gifsicle "gifsicle --version"

echo
echo "... checking additional routines (not urgently needed)"

CheckCommand latex "latex -v"
CheckCommand dvips "dvips -v"
CheckCommand python "python -V"

CheckCommand ssh "ssh -V"
CheckCommand at "at -l"
CheckCommand dnotify "dnotify --version"

CheckCommand "Acrobat Reader" "acroread -help"

if [ -n "$missing" ]; then
  echo ""
  echo "The following programs seem not be installed: "
  echo ""
  echo    "$missing"
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
fi

# xanim
