#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

my $dx = 1000.;
my $dy = 1000.;

$grid->readgrd($file);

translate($grid,$dx,$dy);

$grid->writegrd("new_transl.grd");

#######################################################################

sub translate
{
  my ($grid,$dx,$dy) = @_;

  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};

    $x += $dx;
    $y += $dy;

    $node->{x} = $x;
    $node->{y} = $y;
  }
  
}

