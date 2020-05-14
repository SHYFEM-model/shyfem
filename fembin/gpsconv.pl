#!/usr/bin/perl -s -i.bak
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# possible options: -portrait -landscape

while(<>) {

  if( /%%Orientation: / ) {
    if( $portrait ) {
      $_ = "%%Orientation: Portrait\n";
    } elsif( $landscape ) {
      $_ = "%%Orientation: Landscape\n";
    } else {
      ;
    }
  }

  print;
}
