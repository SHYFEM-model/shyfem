#!/usr/bin/perl
#
# changes .f90 files back to .f format (only continuation lines)
#
#------------------------------------------

while(<>) {

  chomp;

  s/\s*\&\s*$//;

  print "$_\n";
}


#------------------------------------------

