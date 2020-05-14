#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# easy script to change ps to eps

$_ = <>;

print "%!PS-Adobe-2.0 EPSF-2.0\n";

while(<>) {

  if( /^%%EndComments/ ) {
    print;
    $_ = "1.000 1.000 scale\n";
  }

  print;

}
  
