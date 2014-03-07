#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $bline = new grd;
$bline->clone_grid($grid);
$bline->grid_info();
$bline->make_bound_line();

$bline->writegrd("bndline.grd");;

#######################################################################

