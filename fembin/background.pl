#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# background mesh creation routine
#
#----------------------------------------------------------------------
#
# version 1.0
#
# 11.04.2012		finished with call to mesh
#
#----------------------------------------------------------------------
#
# Usage: backb.pl file(s)
#
#       file(s) contain lines from which elements are created
#
#----------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use strict;

use grd;
use grdline;
use mtree;

$::lprefix = "ltmp";
$::gprefix = "gtmp";

my $grid = new grd;
my $infile = $ARGV[0];

my $tree = new mtree(\&is_line_in_line,\&get_line_id);

$grid->readgrd($infile);

my $lines = $grid->get_lines();

foreach my $litem (values %$lines) {
  my ($x,$y) = $grid->make_xy($litem);
  my $grdl = new grdline;
  $grdl->set_line($x,$y);
  $litem->{grdl} = $grdl;

  $tree->insert($litem);
  #$tree->print();

  write_file_from_line($grid,$litem);
}

print "starting printing line info...\n";
$tree->print();

$tree->apply(\&mesh_back);

##############################################################

sub is_line_in_line {

  my ($litem1,$litem2) = @_;

  my $n1 = $litem1->{number};
  my $n2 = $litem2->{number};

  #print "checking if $n2 is in $n1\n";

  my $grdl1 = $litem1->{grdl};
  my $grdl2 = $litem2->{grdl};

  my $x = $grdl2->{x};
  my $y = $grdl2->{y};
  my $x0 = $x->[0];
  my $y0 = $y->[0];

  #print "using point $x0 $y0\n";

  my $result = $grdl1->in_line($x0,$y0);

  #print "result: $result\n";

  return $result;
}

sub get_line_id {

  my $litem = shift;

  return $litem->{number};
}

##############################################################

sub mesh_back {

  my ($node,$parent,$children) = @_;

  my $n = $node->{number};

  my $final = "$::gprefix$n.grd";

  my $meshline = "mesh $::lprefix$n";
  my $line = "";
  foreach my $child (@$children) {
    my $nc = $child->{number};
    $meshline .= " $::lprefix$nc";
    $line .= " $nc";
  }

  print "function called for $n with children $line\n";
  print "    (function call: $meshline)\n";
  system("$meshline; mv final.grd $final");
}

sub write_file_from_line {

  my ($grid,$item) = @_;

  my $numb = $item->{number};
  $item->{type} = 2;
  change_depth($grid,$item);

  my $newgrid = new grd;
  $newgrid->clone_line($item);
  $newgrid->clone_needed_nodes($grid);
  change_node_type($newgrid,0);

  my $file = "$::lprefix$numb.grd";
  $newgrid->writegrd($file);
}

sub invert_line {

  my ($item) = @_;

  my $vert = $item->{vert};
  my @new = reverse(@$vert);
  $item->{vert} = \@new;
}

sub change_node_type {

  my ($grid,$type) = @_;

  my $nodes = $grid->get_nodes();

  foreach my $item (values %$nodes) {
    $item->{type} = $type;
  }
}

sub change_depth {

  my ($grid,$item) = @_;

  my $vert = $item->{vert};
  my $depth = $item->{h};

  foreach my $node (@$vert) {
    my $nitem = $grid->get_node($node);
    $nitem->{h} = $depth;
  }

}

##############################################################
