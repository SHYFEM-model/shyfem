#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# substitutes values in given column in file (default: col=0)
#
# -col=#	column to substitute
# -row=#	if given substitute only in this row, else all
# -value=#	value to substitute
#
# columns are counted from 0, rows from 1
#
#-------------------------------------------------------------

use strict;

$::col = 0 unless $::col;
$::row = 0 unless $::row;
$::value = 0 unless $::value;

my $i = 0;

while(<>) {

  chomp;
  s/^\s+//;
  s/,/ /g;
  my @f = split;
  $i++;

  if( $::row == 0 or $::row == $i ) {
    $f[$::col] = $::value;
  }

  $_ = join(" ",@f);

  print "$_\n";
}


