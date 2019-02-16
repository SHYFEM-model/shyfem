#!/usr/bin/perl 

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

$oldpattern = quotemeta('FEMDIR=$HOME/fem');
$newpattern = 'FEMDIR=${SHYFEMDIR:-$HOME/shyfem}' . "\n";

while(<>) {
  if( /^$oldpattern/ ) {
    $_ = $newpattern;
  }
  print;
}

