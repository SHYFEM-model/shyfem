#!/bin/sh

comp=ifort
dmain=dummy_main.f

#--------------------------------------

MakeDummy()
{
  echo "	end" > $dmain
}

#--------------------------------------

MakeDummy

$comp $dmain $*

rm -f $dmain

