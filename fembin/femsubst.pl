#!/usr/bin/perl 

$oldpattern = quotemeta('FEMDIR=$HOME/fem');
$newpattern = 'FEMDIR=${FEMDIR:-$HOME/fem}' . "\n";

while(<>) {
  if( /^$oldpattern/ ) {
    $_ = $newpattern;
  }
  print;
}

