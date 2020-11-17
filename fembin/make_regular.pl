#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# makes regular grd given parameters on command line
#
#------------------------------------------------

use strict;

my $ntype = 5;
my $ltype = 5;

my $nx = shift;
my $ny = shift;
my $x0 = shift;
my $y0 = shift;
my $dx = shift;
my $dy = shift;

unless( $dy ) {
  die "Usage: make_regular.pl nx ny x0 y0 dx dy\n";
}

my $n = 0;
my $l = 0;

#-------------------------------------------------
# write nodes
#-------------------------------------------------

for( my $iy=0; $iy<$ny; $iy++ ) {
  my $y = $y0 + $iy * $dy;
  for( my $ix=0; $ix<$nx; $ix++ ) {
    my $x = $x0 + $ix * $dx;
    $n++;
    print "1 $n $ntype $x $y\n";
  }
}

#-------------------------------------------------
# write lines
#-------------------------------------------------

for( my $ix=0; $ix<$nx; $ix++ ) {
  my $istart = $ix + 1;
  $l++;
  print "3 $l $ltype $ny";
  for( my $iy=0; $iy<$ny; $iy++ ) {
    print "\n" unless $iy%5;
    my $ip = $istart + $iy*$nx;
    print " $ip";
  }
  print "\n";
}
  
for( my $iy=0; $iy<$ny; $iy++ ) {
  my $istart = $nx*$iy + 1;
  $l++;
  print "3 $l $ltype $nx";
  for( my $ix=0; $ix<$nx; $ix++ ) {
    print "\n" unless $ix%5;
    my $ip = $istart + $ix;
    print " $ip";
  }
  print "\n";
}
  
#-------------------------------------------------
# end of routine
#-------------------------------------------------

