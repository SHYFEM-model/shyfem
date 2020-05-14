#!/usr/bin/perl -w -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use grdline;
use strict;

#---------------------------------------------------- options -----------

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

#-------------------------------------------------- main ----------------

my $lines = $grid->get_lines();
my $flag = {};

foreach my $line (values %$lines) {
  my $grline = make_grdline($grid,$line);	#sets up new datastructure
}

$grid->writegrd("smooth.grd");

#-----------------------------------------------------------------

sub make_grdline {

  my ($grid,$line) = @_;

  my ($x,$y) = $grid->make_xy($line);
  my $grline = new grdline;
  $grline->set_line($x,$y);

  return $grline;
}

#-----------------------------------------------------------------

