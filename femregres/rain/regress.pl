#!/usr/bin/perl

$thresh = 1.e-5;

while(<>) {

  chomp;
  s/^\s+//;

  ($time,$val) = split;

  if( $time == 43200 ) {
    $tot += abs($val-0.05);
    $n++;
  } elsif( $time == 86400 ) {
    $tot += abs($val-0.10);
    $n++;
  }
}

$med = $tot / $n if $n > 0;

print "$n  $tot  $med\n\n";

if( $med < $thresh ) {
  print "test passed\n";
} else {
  die "*** test not passed ... $med  $thresh\n";
}
