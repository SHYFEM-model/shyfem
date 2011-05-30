#!/usr/bin/perl -s

$_ = <>;
chomp;

@f = split;
$n = @f;

if( $noextra ) {
  s/\w+\s*$// if $n >= 5;
} elsif( $version ) {
  $_ = $f[1];
} elsif( $date ) {
  $_ = $f[2];
} elsif( $tag ) {
  $_ = $f[3];
} elsif( $extra ) {
  $_ = $f[4];
}

print "$_\n";

