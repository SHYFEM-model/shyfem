#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

while(<>) {

  if( /^%%Pages:\s+(\d+)/ ) {
    $pages = $1;
    die "$pages\n";
  }
}

print "0\n";

