#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# gets name of subroutines/functions defined in file


while(<>) {

	chomp;

	if( /^\s+subroutine\s+(\w+)/i ) {
	  print "$1\n";
	} elsif( /^\s+function\s+(\w+)/i ) {
	  print "$1\n";
	} elsif( /^\s+([\w\s]+)function\s+(\w+)/i ) {
	  $name = $2;
	  $mod = $1;
	  $mod =~ s/\s+//g;
	  if( $mod =~ /^(integer|real|logical|character|doubleprecision)$/i ) {
		print "$name\n";
	  }
	}
}

