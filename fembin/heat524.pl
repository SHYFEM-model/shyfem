#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# rewrites old heat file (5 columns) to new one (4 columns)
#
#--------------------------------------------------------

while(<>) {

  chomp;
  s/^\s+//;

  @f = split;

  $f[4] = $f[5];
  pop(@f);

  $line = join(" ",@f);

  print "$line\n";
}

