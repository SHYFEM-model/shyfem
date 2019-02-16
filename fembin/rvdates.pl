#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# gets revision dates from fortran file

$date = "1998";

if( @ARGV[0] eq "-d" ) {
  $date = @ARGV[1];
  shift; shift;
}

while(<>) {

  if( /^c \s*\$Id/ ) {
	;
  } elsif( /^c .*$date/ ) {
	s/\s+/ /g;
	s/^c //;
	print "$ARGV: $_\n";
  }
}

