#!/usr/bin/perl -ws
#
# translates a grid file
#
# possible command line arguments:
#
# -x0=#  	new origin in x direction (default 0)
# -y0=#  	new origin in y direction (default 0)
# -scale=#	scale with this number (default 1)
# -center	use center of mass of points for translation
# -proj		make projection from lat/lon to cart or viceversa (using x0,y0)
#		(x0,y0 must be lon/lat of center of projection)
#
#----------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

$::scale = 1. unless $::scale;
$::x0 = 0. unless $::x0;
$::y0 = 0. unless $::y0;
$::center = 0 unless $::center;
$::proj = 0 unless $::proj;

#-------------------------------------------------

my $grid = new grd;
my $file = $ARGV[0];

my $x0 = $::x0;
my $y0 = $::y0;
my $scale = $::scale;

$grid->readgrd($file);

if( $::proj ) {
  project($grid,$x0,$y0);
} else {
  ($x0,$y0) = get_center($grid) if $::center;
  translate($grid,$scale,$x0,$y0);
}

$grid->writegrd("new_transl.grd");

#--------------------------------------------------------

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

