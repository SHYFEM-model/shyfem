#!/usr/bin/perl 

$oldpattern = quotemeta('FEMDIR=$HOME/fem');
$newpattern = 'FEMDIR=${SHYFEMDIR:-$HOME/shyfem}' . "\n";

while(<>) {
  if( /^$oldpattern/ ) {
    $_ = $newpattern;
  }
  print;
}

