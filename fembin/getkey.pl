#!/usr/bin/perl -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

if( $h or $ $help or not $ARGV[0] ) {
  die "Usage: getkey.pl keyword file\n";
}

$keyword = shift;
$key = quotemeta($keyword);

#print STDERR "keyword: $keyword - $key\n";

while(<>) {

  if( s/^\s*$key\s*:?\s*// ) {
    print;
  }

}

