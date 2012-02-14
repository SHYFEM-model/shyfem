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
#
#----------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

$::scale = 1. unless $::scale;
$::x0 = 0. unless $::x0;
$::y0 = 0. unless $::y0;
$::center = 0 unless $::center;

#-------------------------------------------------

my $grid = new grd;
my $file = $ARGV[0];

my $x0 = $::x0;
my $y0 = $::y0;
my $scale = $::scale;

$grid->readgrd($file);

($x0,$y0) = get_center($grid) if $::center;
translate($grid,$scale,$x0,$y0);

$grid->writegrd("new_transl.grd");

#--------------------------------------------------------

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

