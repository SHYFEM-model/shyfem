#!/usr/bin/perl -w -s

use lib "$ENV{HOME}/fem/femlib/perl";
use grd;
use strict;

#---------------------------------------------------- options -----------

while( my $file = shift @ARGV ) {

my $grid = new grd;
#my $file = $ARGV[0];
$grid->{verbose} = 0;

$grid->readgrd($file);				#FEM grid

my $elems = $grid->get_elems();
my $ne = scalar values %$elems;
my $nodes = $grid->get_nodes();
my $nn = scalar values %$nodes;
my $lines = $grid->get_lines();
my $nl = scalar values %$lines;

print STDERR "$file:   nodes $nn    elems $ne    lines $nl\n";

}

