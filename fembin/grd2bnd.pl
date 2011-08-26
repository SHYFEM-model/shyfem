#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $lines = $grid->get_lines();

foreach my $line (values %$lines) {
  my $flag = 1;
  my $vert = $line->{vert};
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $x = $node->{x};
    my $y = $node->{y};
    print "$x $y   $flag\n";
    $flag = 0;
  }
}

