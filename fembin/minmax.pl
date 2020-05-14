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
# computes min/max of given column in file (default: col=0)
#
# -col=col	# use column col, else use 0
# -diag		# output can be used to plot diagonal for scatterplot
# -min|-max	# only print min/max
#
#-------------------------------------------------------------

use strict;

$::col = 0 unless $::col;
$::min = 0 unless $::min;
$::max = 0 unless $::max;
$::diag = 0 unless $::diag;
$::h = 0 unless $::h;
$::help = 0 unless $::help;

my $nfiles = @ARGV;
if( $::h or $::help or $nfiles == 0 ) {
  Usage(); exit 0;
}

my $file = $ARGV[0];
if ( not -e $file ) {
  die "File is not existing: $file\n";
}

my ($min,$max);

while(<>) {

  chomp;
  s/^\s+//;
  s/,/ /g;
  my @f = split;

  my $val = $f[$::col];

  $min = $val unless defined $min;
  $min = $val if $val < $min;
  $max = $val unless defined $max;
  $max = $val if $val > $max;

}

if( $::diag ) {
  print "$min $min\n";
  print "$max $max\n";
} elsif( $::min ) {
  print "$min\n";
} elsif( $::max ) {
  print "$max\n";
} else {
  print "$min $max\n";
}

#------------------------------------------------------

sub Usage
{
  print "Usage: minmax.pl [-h|-help] [-options] file(s)\n";
  print "  prints min/max of specified column of file\n";
  print "  options:\n";
  print "  -h|-help         this help screen\n";
  print "  -col=col         use column col to compute min/max (default 0)\n";
  print "  -min             only print minimum\n";
  print "  -max             only print maximum\n";
  print "  -diag            print diagonal to be used in scatterplot\n";
}

#------------------------------------------------------

