#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# translates a grid file
#
# possible command line arguments:
#
# -x0=#  	new origin in x direction (default 0)
# -y0=#  	new origin in y direction (default 0)
# -scale=#	scale with this number (default 1)
# -angle=#	angle for rotation around (x0,y0) or center in degrees
# -center	use center of mass of points for translation
# -proj		make projection from lat/lon to cart or viceversa (using x0,y0)
#		(x0,y0 must be lon/lat of center of projection)
# -no_sort	preserves order of node/elem numbering
#
#----------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;
$::help = 1 if $::help;
$::h = 1 if $::h;
$::scale = 1. unless $::scale;
$::angle = 0. unless $::angle;
$::x0 = 0. unless $::x0;
$::y0 = 0. unless $::y0;
$::center = 0 unless $::center;
$::proj = 0 unless $::proj;
$::no_sort = 0 unless $::no_sort;

if( $::h or $::help ) {
  FullUsage();
    exit 0;
} elsif( not $ARGV[0] ) {
  Usage();
    exit 1;
  print STDERR "No line given -> applying selection to all items\n";
}

#-------------------------------------------------

my $grid = new grd;
my $file = $ARGV[0];

my $x0 = $::x0;
my $y0 = $::y0;
my $scale = $::scale;
my $angle = $::angle;

$grid->readgrd($file);

$grid->set_preserve_order(1) if $::no_sort;

if( $::proj ) {
  project($grid,$x0,$y0);
} else {
  ($x0,$y0) = get_center($grid) if $::center;
  if ($::angle){
    rotate($grid,$angle,$x0,$y0);
  } else {
    translate($grid,$scale,$x0,$y0);
  }
}

$grid->writegrd("new_transl.grd");

#------------------------------------------------------------------------

sub FullUsage {
  print STDERR "                                    \n";
  Usage();
  print STDERR "                                    \n";
  print STDERR "  translates, scales or rotates a grid file\n";
  print STDERR "                                    \n";
  print STDERR "  -h|-help      this help screen\n";
  print STDERR "                                    \n";
  print STDERR " -x0=#  	new origin or center of rotation in x (default 0)\n";
  print STDERR " -y0=#  	new origin or center of rotation in y (default 0)\n";
  print STDERR " -scale=#	scale factor (default 1)\n";
  print STDERR " -angle=#	angle for rotation around (x0,y0) in degrees (default 0)\n";
  print STDERR " -center	use center of mass for translation or rotation\n";
  print STDERR " -proj		make projection from lat/lon to cart or viceversa (using x0,y0)\n";
  print STDERR " 		(x0,y0 must be lon/lat of center of projection)\n";
  print STDERR " -no_sort	preserves order of node/elem numbers\n";
}

sub Usage {
  print STDERR "Usage: grd_transl.pl [-h|-help] [-options] grid\n";
}

#----------------------------------------------------------

sub project
{
  my ($grid,$x0,$y0) = @_;

  my $nodes = $grid->get_nodes();

  my $latlon;
  my $pi = atan2(1.,1.);
  my $rad = $pi / 180.;

  my $yfact = 60 * 1852;               #distance of 1 min = 1 nautical mile
  my $xfact = $yfact * cos($y0*$rad);

  print STDERR "using x0 = $x0 and y0 = $y0\n";
  if( is_latlon($nodes) ) {
    $latlon = 1;
    print STDERR "projecting from latlon to cartesian...\n";
  } else {
    $latlon = 0;
    print STDERR "projecting from cartesian to latlon ...\n";
  }

  foreach my $node (values %$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};

    if( $latlon ) {
      $x = ($x-$x0)*$xfact;
      $y = ($y-$y0)*$yfact;
    } else {
      $x = $x0 + ($x/$xfact);
      $y = $y0 + ($y/$yfact);
    }

    $node->{x} = $x;
    $node->{y} = $y;
  }
  
}

sub translate
{
  my ($grid,$scale,$x0,$y0) = @_;

  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};

    $x -= $x0;
    $y -= $y0;

    $node->{x} = $scale * $x;
    $node->{y} = $scale * $y;
  }
  
}

sub rotate 
{
  my ($grid,$angle,$x0,$y0) = @_;
  
  my $rangle = $angle*3.141592653589/180;	
  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $x = $node->{x} - $x0;
    my $y = $node->{y} - $y0;

    my $x1 = $x*cos($rangle) - $y*sin($rangle);
    my $y1 = $x*sin($rangle) + $y*cos($rangle);
    
    $node->{x} = $x1 + $x0;
    $node->{y} = $y1 + $y0;
  }

}

#--------------------------------------------------------

sub get_center {

  my ($grid) = @_;

  my ($xm,$ym,$n) = 0;

  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};

    $n++;
    $xm += $x;
    $ym += $y;
  }
  
  if( $n ) {
    $xm /= $n;
    $ym /= $n;
  }

  return ($xm,$ym);
}

#--------------------------------------------------------

sub is_latlon {

  my $nodes = shift;

  foreach my $node (values %$nodes) {
    my $x = abs($node->{x});
    my $y = abs($node->{y});

    if( $x > 360 or $y > 90 ) {
      return 0;
    }
  }
  return 1;
}

#--------------------------------------------------------

