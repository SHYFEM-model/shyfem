#!/usr/bin/perl -w
#
# writes node list to GRD format
#
#----------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];
my $nodes = $ARGV[1];
my @nodes = ();
my $n = 0;

$grid->readgrd($file);

open(FILE,"<$nodes");

while(<FILE>) {

  chomp;
  my $nnode = $_;
  last if $nnode == 0;

  my $node = $grid->get_node($nnode);
  my $x = $node->{x};
  my $y = $node->{y};
  print "1 $nnode 0 $x $y\n";

  push(@nodes,$nnode);
}

my $nn = @nodes;

print "\n";
print "3 1 0 $nn\n";

my $i = 0;
foreach my $nnode (@nodes) {
  $i++;
  print " $nnode";
  print "\n" if $i%8 == 0;
}
print "\n";

