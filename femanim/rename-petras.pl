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
# Script changes the file names like *.n.* (n are numbers 1,2,3..)   
# to names, where n has the same number of digits, by adding leading zeros.
#
#---------------------------------------------------------------------

$num_frame = scalar @ARGV;

$zero = '0000000000000000';

$ndig = length($num_frame);

while( $file = shift )
{
  @f = split(/\./,$file);
  $n = $f[1];  
  $nzero = $ndig - length($n);  
  $zzero = substr($zero, 1, $nzero);  
  $num = $zzero . $n;

  $f[1] = $num;
  $newfile = join(".",@f);  
  print "$file -> $newfile\n";
  rename $file , $newfile;
}

