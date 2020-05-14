#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# finds subroutines and functions in files
#
#---------------------------------------------------

while(<>) {

  chomp;

  if( /^[cC*]/ or /\s*!/ ) { next; }

  if( /\s+(subroutine)\s+(\w+)/i or /\s+(function)\s+(\w+)/i ) {
    print "$1   $2\n";
  }
}

#---------------------------------------------------

