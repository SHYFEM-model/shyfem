#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

while(<>) {

  chomp;

  $a = ord('a');
  $z = ord('z');
  $A = ord('A');
  $Z = ord('Z');

  $len = length($_);

  for( $l=0; $l<$len; $l++ ) {
    $c = substr($_,$l,1);
    $ic = ord($c);
    #print "$c $ic $a $z\n";
    if( $a <= $ic and $ic <= $z ) {
      $il++;
    } elsif( $A <= $ic and $ic <= $Z ) {
      $iu++;
    } else {
      $io++;
    }
    $it++;
  }
}

$if = $iu / ($iu+$il);
#print "if=$if iu=$iu il=$il io=$io it=$it\n";
print "$if\n";

if( $if < 0.1 ) {
  exit 0;
} else {
  exit 1;
}
