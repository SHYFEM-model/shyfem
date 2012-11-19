#!/usr/bin/perl -ws
#
# computes min/max of given column in file (default: col=0)
#
# -col=#
# -diag		# output can be used to plot diagonal for scatterplot
#
#-------------------------------------------------------------

use strict;

$::col = 0 unless $::col;
$::diag = 0 unless $::diag;

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
} else {
  print "$min $max\n";
}

