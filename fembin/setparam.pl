#!/usr/bin/perl 

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

####!/usr/bin/perl -i.bak

$param = shift;
$value = shift;

while(<>) {
  if( /^\s+parameter/ ) {
    if( /(\W$param\s*=\s*)([+-]?\d+)/ ) {
        print STDERR "$1 -> $2\n";
    }
  }
}

