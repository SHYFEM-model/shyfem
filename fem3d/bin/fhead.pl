#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# gets function headers and comments

while(<>) {

  chomp;

  if( $in_proc ) {
    #if( /^[c!*]\s*(\S.*)/ ) {	# comment
    if( /^c\s*(\S.*)/ ) {	# comment
      $comm = $1;
      $in_proc = 0;
	#print "with comment... $_\n";
	#print "with comment... $1\n";
    } elsif( /^\s+\S/ ) {	# no comment
      $comm = "";
      $in_proc = 0;
	#print "without comment... $_\n";
    }
    unless( $in_proc ) {
      print "c $first\n";
      print "c\t$comm\n" if $comm;
    }
  }

  if( /^\s+function/ || /^\s+subroutine/ ) {
    s/^\s+//;
    $first = $_;
    $in_proc = 1;
  }

}
