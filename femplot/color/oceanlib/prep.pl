#!/usr/bin/perl

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
