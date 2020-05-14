#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# changes macros in Makefile
#
#----------------------------------------------

while(<>) {

  chomp;

  change("FORTRAN_COMPILER",$FORTRAN_COMPILER);
  change("C_COMPILER",$C_COMPILER);
  change("PARALLEL_OMP",$PARALLEL_OMP);
  change("SOLVER",$SOLVER);
  change("NETCDF",$NETCDF);
  change("NETCDFDIR",$NETCDFDIR);
  change("GOTM",$GOTM);
  change("ECOLOGICAL",$ECOLOGICAL);

  print "$_\n";
}

sub change
{
  my ($name,$new_value) = @_;

  return unless $new_value;

  if( /^\s*$name\s*=\s*(\S+)/ ) {
    #print STDERR "$name | $new_value | $1\n";
    return if $1 eq $new_value;
    $_ = "#$_\n" . "$name = $new_value";
  }
}

