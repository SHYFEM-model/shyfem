#!/usr/bin/perl
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

