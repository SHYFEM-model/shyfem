#!/usr/bin/perl

# gets keywords from header of fortran file

while(<>) {

  if( /^[cC\!]\s+(.*)\s+:\s*\n$/ ) {
    print "$1\n";
  }

}
