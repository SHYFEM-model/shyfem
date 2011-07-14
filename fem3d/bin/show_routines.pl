#!/usr/bin/perl -n

  $no_space = $_;
  $no_space =~ s/\s+//g;

  print if $no_space =~ /^subroutine/i;
  print if $no_space =~ /^function/i;

  print if $no_space =~ /^(real|integer|doubleprecision|logical)function/i;

