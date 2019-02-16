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
# prints out function and subroutine definitions
#
#-----------------------------------------------

while(<>) {

  if( /function/ ) {
    check_function($1,$_);
  } elsif( /subroutine/ ) {
    check_function($1,$_);
  }
}

sub check_function {

  my ($key,$line) = @_;

  return if( $line =~ /^\s+end/ );
  return if( $line =~ /^[cC*]/ );
  return if( $line =~ /^\s*\!/ );

  print_to_end($line);
}

sub print_to_end {

  my ($line) = @_;

  print $line;
  return unless /\(/;			# if no "(" we are finished

  until( $line =~ /\)\s*$/ ) {		# ")" must be last char on line
    $line = <>;
    print $line;
  }
}

