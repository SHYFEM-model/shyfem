#!/usr/bin/perl -w

my $n = 0;

while(<>) {

  chomp;
  $n++ if( /VERS_/ );

  last if $n > 1;

  print "$_\n";
}
