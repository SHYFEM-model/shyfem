#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# gets keywords from header of fortran file

while(<>) {

  if( /^[cC\!]\s+(.*)\s+:\s*\n$/ ) {
    print "$1\n";
  }

}
