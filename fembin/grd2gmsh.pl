#!/usr/bin/perl -w
#
# convert line in GRD format into GMSH format
#
#-------------------------------------------------------------

use lib "$ENV{HOME}/lib/perl";
use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $lines = $grid->get_lines();

foreach my $line (values %$lines) {
  my $nvert = $line->{nvert};
  my $flag = 1;
  print "$nvert\t$flag\n";
  my $vert = $line->{vert};
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $x = $node->{x};
    my $y = $node->{y};
    print "$x\t$y\n";
  }
}

