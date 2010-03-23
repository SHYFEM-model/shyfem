#!/usr/bin/perl -s -i.bak
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
