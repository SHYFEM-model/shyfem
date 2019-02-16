#!/usr/bin/perl -s
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# modifies parameters in files created with gp
#
# modifies file inline
#
# options:
#		-h|-help
#		-bw   -color
#		-width=#
#
#--------------------------------------------------------

if( $h or $help ) {
  print STDERR "Usage: gp_modify.pl [-h|-help] [-options] file\n";
  print STDERR "   options:\n";
  print STDERR "     -h|-help     this help screen\n";
  print STDERR "     -bw          make plot black and white\n";
  print STDERR "     -color       make plot color\n";
  print STDERR "     -width=#     set line width to #\n";
  exit(0);
}

while(<>) {

  chomp;

  if( $width and /^\/gnulinewidth/ ) {
    @f = split;
    $f[1] = $width;
    $_ = join(" ",@f);
  } elsif( $bw and /^\/Color/ ) {
    @f = split;
    $f[1] = "false";
    $_ = join(" ",@f);
  } elsif( $color and /^\/Color/ ) {
    @f = split;
    $f[1] = "true";
    $_ = join(" ",@f);
  }

  print "$_\n";
}
