#!/usr/bin/perl -s -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# finds closest grid nodes to sparse nodes (in grd format)
#
#--------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

$::nosort = 0 unless $::nosort;		# do not pre-sort points
$::silent = 0 unless $::silent;		# try to be silent (no messages)
$::nodeonly = 0 unless $::nodeonly;	# only write node to stdout

#-------------------------------------------------------------

my $grid_file = $ARGV[0];
my $node_list = $ARGV[1];

unless( $node_list ) {
  die "Usage: grd_closest.pl grd-file grd-nodes\n";
}

my $grid = new grd;
$grid->set_verbose(0) if $::silent;
$grid->readgrd("$grid_file");
my $ngrid = new grd;
$ngrid->set_verbose(0) if $::silent;
$ngrid->readgrd("$node_list");

my @list = ();
open(GRD,">closest_nodes.grd");
open(TXT,">closest_nodes.txt");

#------------------------------------------------------------

my @nodes = $ngrid->get_nodes_ordered();
@nodes = sort by_number @nodes unless $::nosort;

foreach my $number ( @nodes ) {
  my $nitem = $ngrid->get_node( $number );
  my $x = $nitem->{x};
  my $y = $nitem->{y};

  my $n = get_closest($grid,$x,$y);
  my $text = "$number    $n";
  $text = "$n" if $::nodeonly;
  push(@list,$text);

  my $item = $grid->get_node($n);
  my $type = $item->{type};
  $number = $item->{number};
  $x = $item->{x};
  $y = $item->{y};
  print GRD "1 $number $type $x $y\n" unless $::nodeonly;
  print TXT "$number\n" unless $::nodeonly;
}

foreach my $line (@list) {
  print "$line\n";
}

unless( $::silent) {
  print STDERR "nodes have been written to closest_nodes.grd";
  print STDERR " and closest_nodes.txt\n";
}

#------------------------------------------------------------

sub get_closest {

  my ($grid,$xn,$yn) = @_;

  my $nn = 0;
  my $dist = 1.e+30;

  my $nodes = $grid->get_nodes();
  foreach my $item (values %$nodes) {
    my $number = $item->{number};
    my $x = $item->{x};
    my $y = $item->{y};

    my $d = ($x-$xn)*($x-$xn) + ($y-$yn)*($y-$yn);
    if( $d < $dist ) {
      $dist = $d;
      $nn = $number;
    }
  }

  return $nn;
}

#------------------------------------------------------------

sub by_node_number {

  if( $a->{number} > $b->{number} ) {
    return 1;
  } elsif( $a->{number} < $b->{number} ) {
    return -1;
  } else {
    return 0;
  }
}

#------------------------------------------------------------

sub by_number {

  if( $a > $b ) {
    return 1;
  } elsif( $a < $b ) {
    return -1;
  } else {
    return 0;
  }
}

#------------------------------------------------------------

