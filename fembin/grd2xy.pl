#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# extracts x/y/z info from grd file and writes it to xy file
#
# z value may be missing
#
# -csv	format of comma seperated value
#
#------------------------------------------------------------------------

while(<>) {

  chomp;
  s/^\s*//;	# get rid of leading spaces

  @f = split;
  $nf = @f;

  next unless $f[0] == 1;	#only transform nodes

  if( $nf == 5 ) {
    $x = $f[3];
    $y = $f[4];
    $z = "";
  } elsif( $nf == 6 ) {
    $x = $f[3];
    $y = $f[4];
    $z = $f[5];
  } else {
    die "expecting 5 or 6 values: $_\n";
  }

  if( $csv ) {
    print "$x;$y;$z\n";
  } else {
    print "$x $y $z\n";
  }
}

#------------------------------------------------------------------------

