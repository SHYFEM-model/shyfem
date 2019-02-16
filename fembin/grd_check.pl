#!/usr/bin/perl -w -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use grdline;
use strict;

#---------------------------------------------------- options -----------

$::help = 1 if $::help;
$::h = 1 if $::h;
$::invert = 0 unless $::invert;
$::exclude = 0 unless $::exclude;

if( $::n or $::e or $::l ) {	#explicitly given -> only on these items
  $::n = 0 unless $::n;
  $::e = 0 unless $::e;
  $::l = 0 unless $::l;
} else {			#default -> on all items
  $::n = 1;
  $::e = 1;
  $::l = 1;
}

$::print = 0 unless $::print;
$::delete = 0 unless $::delete;
$::type = -1 unless defined($::type);

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
}

print STDERR "in: $::inside  out: $::outside  exclude: $::exclude\n";
print STDERR "print: $::print  type: $::type  delete: $::delete\n";
print STDERR "n: $::n  e: $::e  l: $::l\n";

#------------------------------------------------- read files -----------

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);				#FEM grid

exit;

#-------------------------------------------------- main ----------------

my $flag = {};
my $lines = {};
foreach my $line (values %$lines) {

  loop_on_elements($grid,$flag,\&modify_element) if $::e;
  loop_on_lines($grid,$flag,\&modify_line) if $::l;
  loop_on_nodes($grid,$flag,\&modify_node) if $::n;
}

$grid->writegrd("modify.grd");

#------------------------------------------------------------------------

sub FullUsage {
  print STDERR "                                    \n";
  Usage();
  print STDERR "                                    \n";
  print STDERR "  Selects and modifies items of grid contained in line\n";
  print STDERR "                                    \n";
  print STDERR "  -h|-help     this help screen\n";
  print STDERR "                                    \n";
  print STDERR "  {-n|-e|-l}   modify nodes,elements,lines (default: all)\n";
  print STDERR "  -print       print selected items\n";
  print STDERR "  -type=type   set type of selected items to type\n";
  print STDERR "  -delete      delete selected items\n";
  print STDERR "                                    \n";
  print STDERR "  -invert      invert selection - modify outer area\n";
  print STDERR "  -exclude     exclude border line items from selection\n";
}

sub Usage {
  print STDERR "Usage: grd_modify.pl [-h|-help] [-options] grid line\n";
}

#----------------------------------------------------------

sub loop_on_nodes {

  my ($grid,$flag,$proc) = @_;

  $grid->make_used();
  my $nodes = $grid->get_nodes();

  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $f = $$flag{$number};
    if( $f and $::inside or not $f and not $::inside ) {
      &$proc($grid,$node);
    }
  }
}

sub loop_on_elements {

  my ($grid,$flag,$proc) = @_;

  my $elems = $grid->get_elems();

  loop_on_items($grid,$flag,$proc,$elems);
}
 
sub loop_on_lines {

  my ($grid,$flag,$proc) = @_;

  my $lines = $grid->get_lines();

  loop_on_items($grid,$flag,$proc,$lines);
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

#    if( $n ) {
#	print "flagged: element $number ($n)\n";
#    }

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

sub modify_node {

  my ($grid,$node) = @_;

  my $number = $node->{number};

  if( $::print ) {
    print "modifying node $number\n";
  } elsif( $::type >= 0 ) {
    $node->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_node($node) unless $node->{used};
  }
}

sub modify_element {

  my ($grid,$elem) = @_;

  my $number = $elem->{number};

  if( $::print ) {
    print "modifying element $number\n";
  } elsif( $::type >= 0 ) {
    $elem->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_elem($elem);
  }
}

sub modify_line {

  my ($grid,$line) = @_;

  my $number = $line->{number};

  if( $::print ) {
    print "modifying line $number\n";
  } elsif( $::type >= 0 ) {
    $line->{type} = $::type;
  } elsif( $::delete ) {
    $grid->delete_line($line);
  }
}

#-----------------------------------------------------------------

sub print_flaged_nodes {

  my ($flag) = @_;

  foreach my $number (keys %$flag) {
    print "$number\n";
  }
}

sub flag_nodes {

  my ($grid,$grline,$flag) = @_;

  my $nodes = $grid->get_nodes();
  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $x = $node->{x};
    my $y = $node->{y};
    if( $grline->in_line($x,$y) ) {
      $$flag{$number} = 1;
    }
  }

  #return $flag;
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

