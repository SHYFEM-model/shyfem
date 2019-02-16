#!/usr/bin/perl -s -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# converts node numbers to grd-nodes
#
# possible command line options:
#
#	-line	make line from nodes
#
#--------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use str;
use grd;
use strict;

#-------------------------------------------------------------
# command line options
#-------------------------------------------------------------
$::line = 1 if $::line;
#-------------------------------------------------------------

#-------------------------------------------------------------
#-------------------------------------------------------------

my $grdfile = $ARGV[0];
my $nodefile = $ARGV[1];

unless ( $nodefile ) {
  die "Usage: nodes2grd.pl grd-file node-file\n";
}

my $grid = new grd;
$grid->readgrd("$grdfile");

my $list = read_nodes("$nodefile");

write_nodes($grid,$list);

#------------------------------------------------------------

sub read_nodes {

  my $file = shift;

  my @list = ();

  open(FILE,"<$file") || die "Cannot open file: $file\n";;

  while(<FILE>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    push(@list,$_) if $_;
  }

  close(FILE);

  return \@list;
}

#------------------------------------------------------------

sub write_nodes {

  my ($grid,$list) = @_;

  my $ibtyp = 6;

  open(GRD,">nodes.grd");

  write_nodes_to_grd ($ibtyp,$grid,$list);

  close(GRD);
  print STDERR "nodes written to file nodes.grd\n";
}

sub write_nodes_to_grd {

  my ($ibtyp,$grid,$list) = @_;

  return if( $ibtyp <= 0 );
  my $type = $ibtyp;

  foreach my $node (@$list) {
    if( $grid->exists_node($node) ) {
      my $item = $grid->get_node($node);
      my $x = $item->{x};
      my $y = $item->{y};
      print GRD "1 $node $type $x $y\n";
    } else {
      print STDERR "*** Node $node is not existing ... not converted to grd\n";
    }
  }

  return unless $::line;

  my $nnodes = @$list;

  print GRD "";
  print GRD "3 1 $type $nnodes\n";

  my $i = 0;
  foreach my $node (@$list) {
    print GRD " $node";
    $i++;
    print GRD "\n" if $i%10 == 0;
  }
  print GRD "\n\n";
}


