#!/usr/bin/perl -n

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

  $no_space = $_;
  $no_space =~ s/\s+//g;

  print if $no_space =~ /^subroutine/i;
  print if $no_space =~ /^function/i;

  print if $no_space =~ /^(real|integer|doubleprecision|logical)function/i;

