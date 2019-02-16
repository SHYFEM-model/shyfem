#!/usr/bin/perl -w -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# version = 2.0
#
# 07.10.2010	ggu	act on all items if no line is given
# 07.10.2010	ggu	translate nodes
#
#----------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use grdline;
use strict;

#---------------------------------------------------- options -----------

$::help = 1 if $::help;
$::h = 1 if $::h;
$::invert = 0 unless $::invert;
$::exclude = 0 unless $::exclude;
$::preserve = 0 unless $::preserve;

$::e = 1 if $::e;	#extract

my $stype = 888;	#special type

$::inside = 1;		#if true delete inside
$::outside = 0;		#if true delete outside

if( $::invert ) {
  $::outside = 1;
  $::inside = 0;
}

if( $::h or $::help ) {
  FullUsage();
  exit 0;
} elsif( not $ARGV[0] ) {
  Usage();
  exit 1;
} elsif( not $ARGV[1] ) {
  #Usage();
  #exit 1;
  print STDERR "No line given -> applying selection to all items\n";
}


#------------------------------------------------- read files -----------

my $grid = new grd;
my $file = $ARGV[0];
$grid->readgrd($file);				#FEM grid

my $gline = new grd;
my $gfile = $ARGV[1];
if( $gfile ) {
  $gline->readgrd($gfile);			#grid with elements
} else {
  make_dummy_line($grid,$gline);		# dummy line around all nodes
}

#-------------------------------------------------- main ----------------


if( $::e ) {	#extract macro elements -> delete everything else

  print "extracting macro elements with type $stype\n";
  delete_item_no_type($grid,$grid->get_elems(),$stype);
  delete_item_no_type($grid,$grid->get_lines(),$stype);
  $grid->delete_unused();

} else {	#change type

  my $elems = $gline->get_elems();

  foreach my $elem (values %$elems) {
    my $grline = make_grdline($gline,$elem);	#sets up new datastructure
    my $type = $elem->{type};
    print "changing type to $type\n";
    change_elem_types($grid,$grline,$type);
  }
}

$grid->set_preserve_order($::preserve);
$grid->writegrd("modify.grd");

#------------------------------------------------------------------------

sub delete_item_no_type {

  my ($grid,$items,$stype) = @_;

  foreach my $item (values %$items) {
    my $verts = $item->{vert};
    my $bgood = 0;
    foreach my $nnode (@$verts) {
      my $node = $grid->get_node($nnode);
      my $type = $node->{type};
      if( $type == $stype ) {
        $bgood = 1;
        last;
      }
    }
    if( not $bgood ) {
      my $n = $item->{number};
      delete $$items{$n};
    }
  }

}

#------------------------------------------------------------------------

sub FullUsage {
  print STDERR "                                    \n";
  Usage();
  print STDERR "                                    \n";
  print STDERR "  Selects and modifies items of grid contained in line\n";
  print STDERR "                                    \n";
  print STDERR "  If no line is given acts on all items of grid\n";
  print STDERR "                                    \n";
  print STDERR "  -h|-help     this help screen\n";
  print STDERR "                                    \n";
  print STDERR "  {-n|-e|-l}   modify nodes,elements,lines (default: all)\n";
  print STDERR "  -print       print selected items\n";
  print STDERR "  -type=type   set type of selected items to type\n";
  print STDERR "  -depth=depth set depth of selected items to depth\n";
  print STDERR "  -trans=dx,dy translate selected nodes by dx/dy\n";
  print STDERR "  -delete      delete selected items\n";
  print STDERR "                                    \n";
  print STDERR "  -invert      invert selection - modify outer area\n";
  print STDERR "  -exclude     exclude border line items from selection\n";
  print STDERR "  -preserve    preserves order of items\n";
}

sub Usage {
  print STDERR "Usage: grd_modify.pl [-h|-help] [-options] grid [line]\n";
}

#----------------------------------------------------------

sub loop_on_elements {

  my ($grid,$flag,$proc) = @_;

  my $elems = $grid->get_elems();

  loop_on_items($grid,$flag,$proc,$elems);
}
 
sub loop_on_items {

  my ($grid,$flag,$proc,$items) = @_;

  my $elems = $grid->get_elems();
  foreach my $item (values %$items) {
    my $vert = $item->{vert};
    my $nvert = $item->{nvert};
    my $number = $item->{number};
    my $n = 0;
    foreach my $node (@$vert) {
      $n++ if $$flag{$node};
    }

    if( $::inside ) {
      if( $::exclude ) {
        if( $n == $nvert ) {
	  &$proc($grid,$item);
	}
      } else {
        if( $n ) {
	  &$proc($grid,$item);
	}
      }
    } else {
      if( $::exclude ) {
        if( $n == 0 ) {
	  &$proc($grid,$item);
	}
      } else {
        if( $n < $nvert ) {
	  &$proc($grid,$item);
	}
      }
    }
  }
}

#-----------------------------------------------------------------

sub change_area {

  my ($grid,$elem) = @_;

  my $number = $elem->{number};

  $elem->{type} = $::type;
  print "changing elem: $number $::type\n";
}

#-----------------------------------------------------------------

sub change_elem_types {

  my ($grid,$grline,$type) = @_;

  my $elems = $grid->get_elems();
  foreach my $elem (values %$elems) {
    my ($xm,$ym) = $grid->make_central_point($elem);
    if( $grline->in_line($xm,$ym) ) {
      $elem->{type} = $type;
    }
  }
}

#-----------------------------------------------------------------

sub make_grdline {

  my ($gline,$line) = @_;

  my ($x,$y) = make_xy_array($gline,$line);
  my $grline = new grdline;
  $grline->set_line($x,$y);

  return $grline;
}

sub make_xy_array {

  my ($grid,$line) = @_;

  my @x = ();
  my @y = ();

  my $vert = $line->{vert};
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $x = $node->{x};
    my $y = $node->{y};
    push(@x,$x);
    push(@y,$y);
  }

  return (\@x,\@y);
}

#-----------------------------------------------------------------

sub make_dummy_line {

  my ($grid,$gline) = @_;

  my ($xmin,$ymin,$xmax,$ymax) = $grid->get_xy_minmax();

  print STDERR "xy_minmax: $xmin $ymin $xmax $ymax \n";

  $xmin -= 0.1*($xmax-$xmin);
  $xmax += 0.1*($xmax-$xmin);
  $ymin -= 0.1*($ymax-$ymin);
  $ymax += 0.1*($ymax-$ymin);

  $gline->make_node(1,0,$xmin,$ymin);
  $gline->make_node(2,0,$xmax,$ymin);
  $gline->make_node(3,0,$xmax,$ymax);
  $gline->make_node(4,0,$xmin,$ymax);

  $gline->make_line(1,0,0,1,2,3,4,1);
}

#-----------------------------------------------------------------

