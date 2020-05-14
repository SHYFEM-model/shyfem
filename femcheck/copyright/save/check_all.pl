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
# cleans output from check_all.sh by deleting directory output
#
#-----------------------------------------------------------

$dir = shift;

while(<>) {

  if( /^$dir/ ) {
    skip_directory();
  }

  print;
}

sub skip_directory
{
  <>;
  while(<>) {
    last if /^-----/;
  }
  exit 0 unless /^-----/;	# EOF read
}

