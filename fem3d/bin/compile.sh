#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

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

