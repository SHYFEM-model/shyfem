#!/usr/bin/perl

$n = 0;
$type = 3;

while(<>) {

  chomp;
  s/^\s*//;
  @f = split;
  $nf = @f;

  if( $nf == 2 ) {
    $x = $f[0];
    $y = $f[1];
    $z = "";
  } elsif( $nf == 3 ) {
    $x = $f[0];
    $y = $f[1];
    $z = $f[2];
  } else {
    die "expecting 2 or 3 values: $_\n";
  }

  $n++;
  print "1 $n $type $x $y $z\n";
}

