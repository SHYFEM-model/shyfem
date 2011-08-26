#!/usr/bin/perl -w -s

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use grdline;
use strict;

$::help = 1 if $::help;
$::h = 1 if $::h;
$::invert = 0 unless $::invert;
$::exclude = 0 unless $::exclude;

$::inside = 0;		#if true delete inside
$::outside = 1;		#if true delete outside

if( $::invert ) {
  $::outside = 0;
  $::inside = 1;
}

if( $::h or $::help ) {
  FullUsage();
  exit 0;
} elsif( not $ARGV[1] ) {
  Usage();
  exit 1;
}

print STDERR "in: $::inside  out: $::outside  ex: $::exclude\n";

my $grid = new grd;
my $file = $ARGV[0];
my $gline = new grd;
my $gfile = $ARGV[1];

$grid->readgrd($file);
$gline->readgrd($gfile);

my $lines = $gline->get_lines();
my $flag = {};

foreach my $line (values %$lines) {
  my ($x,$y) = make_xy_array($gline,$line);
  my $grline = new grdline;
  $grline->set_line($x,$y);
  flag_nodes($grid,$grline,$flag);
  #print_flaged_nodes($flag);
  delete_elements($grid,$flag);
}

$grid->writegrd("delete.grd");

#----------------------------------------------------------

sub FullUsage {
  Usage();
  print STDERR "  Extracts elements contained in line and deletes others.\n";
  print STDERR "  -h|-help     this help screen\n";
  print STDERR "  -invert      invert selection - delete inner area\n";
  print STDERR "  -exclude     exclude border line elements from selection\n";
}

sub Usage {
  print STDERR "Usage: extractgrd.pl [-h|-help] [-invert] [-exclude] ";
  print STDERR "grid line\n";
}

#----------------------------------------------------------

sub delete_elements {

  my ($grid,$flag) = @_;

  my $elems = $grid->get_elems();
  foreach my $elem (values %$elems) {
    my $vert = $elem->{vert};
    my $nvert = $elem->{nvert};
    my $n = 0;
    foreach my $node (@$vert) {
      $n++ if $$flag{$node};
    }

    if( $::inside ) {
      if( $::exclude ) {
        if( $n ) {
          $grid->delete_elem($elem);
	}
      } else {
        if( $n == $nvert ) {
          $grid->delete_elem($elem);
	}
      }
    } else {
      if( $::exclude ) {
        if( $n != $nvert ) {
          $grid->delete_elem($elem);
	}
      } else {
        if( $n == 0 ) {
          $grid->delete_elem($elem);
	}
      }
    }

  }

  delete_flagged_nodes($grid,$flag);
}

sub delete_flagged_nodes {

  my ($grid,$flag) = @_;

  $grid->make_used();

  my $nodes = $grid->get_nodes();
  foreach my $node (values %$nodes) {
    my $number = $node->{number};
    my $f = $$flag{$number};
    if( $f and $::inside or not $f and not $::inside ) {
      unless( $node->{used} ) {
        $grid->delete_node($node);
      }
    }
  }
}

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

