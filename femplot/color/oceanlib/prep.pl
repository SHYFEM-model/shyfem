#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

while(<>) {

  if( /\[\[/ ) {
    $incm = 1;
    $file = $ARGV;
    $name = $file;
    $name =~ s/\.py$//;
    s/cm_data/_${name}_data/;
  }

  if( $incm ) {
    print;
  }

  if( /\]\]/ ) {
    $incm = 0;
  }
}
