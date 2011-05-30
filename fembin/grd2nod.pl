#!/usr/bin/perl -w
#
# writes nodes of lines in GRD to stdout
#
#----------------------------------------------------

use lib "$ENV{HOME}/fem/femlib/perl";
use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $lines = $grid->get_lines();

foreach my $line (values %$lines) {
  my $flag = 1;
  my $vert = $line->{vert};
  my $i = 0;
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $n = $node->{number};
    print " $n";
    $i++;
    print "\n" if $i%10 == 0;
  }
  print "\n";
}

