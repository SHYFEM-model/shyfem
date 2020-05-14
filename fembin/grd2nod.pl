#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# writes nodes of lines in GRD to stdout
#
# options:
#
#	-zero	postpend 0
#	-list	write one node per line
#
#----------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$::h = 1 if $::h;
$::help = 1 if $::h;
$::help = 1 if $::help;
$::zero = 1 if $::zero;
$::list = 1 if $::list;

FullUsage() if $::help;
Usage() unless $file;

$grid->readgrd($file);

my $lines = $grid->get_lines();
my @keys = keys %$lines;
@keys = sort {$a<=>$b} @keys;

foreach my $nline (@keys) {
  my $line = $grid->get_line($nline);
  my $flag = 1;
  my $vert = $line->{vert};
  print STDERR "writing line $nline\n";
  $::numbers = 0;
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $n = $node->{number};
    output($n);
  }
  output(0) if $::zero;

  print "\n" if $::numbers != 0 and not $::list;
  print "\n";
}

sub output {

  my $n = shift;

  my $wrap = 10;

  $::numbers++;
  $::numbers = $::numbers % $wrap;

  if( $::list ) {
    print "$n\n";
  } else {
    print " $n";
    print "\n" if $::numbers == 0;
  }
}

sub Usage {

  print "grd2nod.pl [-h|-help] [options] grd-file\n";
  exit 0;
}
sub FullUsage {

  print "grd2nod.pl [-h|-help] [options] grd-file\n";
  print "  -h|-help         this message\n";
  print "  -zero            postpend 0\n";
  print "  -list            write one node per line\n";
  exit 0;
}


