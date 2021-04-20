#!/usr/bin/perl -w -s
# 
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# converst dm2 files to grd files

#------------------------------------------------------------------------

use strict;

while( <> ) {

  chomp;

  if( /^E3T/ ) {
    my @f = split;
    my $ie = $f[1];
    my $k1 = $f[2];
    my $k2 = $f[3];
    my $k3 = $f[4];
    print "2 $ie 0 3 $k1 $k2 $k3\n";
  } elsif( /^ND/ ) {
    my @f = split;
    my $k = $f[1];
    my $x = $f[2];
    my $y = $f[3];
    my $z = $f[4];
    print "1 $k 0 $x $y $z\n";
  }
}

