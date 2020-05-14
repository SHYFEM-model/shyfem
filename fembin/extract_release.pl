#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

my $n = 0;

while(<>) {

  chomp;
  $n++ if( /VERS_/ );

  last if $n > 1;

  print "$_\n";
}
