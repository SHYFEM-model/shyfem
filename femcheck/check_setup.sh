#!/bin/sh

log=CHECKLOG
rm -f $log
echo "quit" > quit.tmp

errors=NO
missing=""

CheckCommand()
{
  name=$1
  command=$2
  error=$3

  if [ -z "$error" ]; then
    error=0
  fi

  ($command) >> $log 2>&1

  status=$?

  if [ $status -eq $error ]; then
    echo "... $name is installed"
  else
    echo "*** $name is not installed"
    missing="$missing $name"
    errors=YES
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

cat > quit.tmp <<EOI
quit
EOI
}

#---------------------------------------------------

CreateInputFiles

CheckCommand make "make -v"
CheckCommand perl "perl -v"
CheckCommand g77 "g77 -v"
CheckCommand gcc "gcc -v"
CheckCommand bash "bash --version"
CheckCommand X11 "gcc -L/usr/X11/lib -L/usr/X11R6/lib -lXt -lX11 test.c"

CheckCommand ssh "ssh -V"

CheckCommand ImageMagic "mogrify -version"
CheckCommand gnuplot "gnuplot quit.tmp"
CheckCommand "ghostview gv" "gv --version"
CheckCommand ghostscript "gs -v"
CheckCommand gifsicle "gifsicle --version"
CheckCommand "Acrobat Reader" "acroread -help"

CheckCommand latex "latex -v"
CheckCommand dvips "dvips -v"
CheckCommand python "python -V"

CheckCommand at "at -l"
CheckCommand dnotify "dnotify --version"

if [ $errors = "YES" ]; then
  echo ""
  echo "The following programs seem not be installed: "
  echo ""
  echo    "$missing"
  echo ""
  echo "They might not be indispensable but you are advised"
  echo "to install these programs from you Linux Distribution CD."
  echo "How exactly this is done depends a lot on the type of"
  echo "distribution you are using. But you can try and bring up"
  echo "a menu Configuration and from there on Install Software."
  echo "Search for the keyword of the missing file and then install"
  echo "the files on your harddisk."
  echo ""
fi

# xanim
