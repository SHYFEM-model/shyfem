#!/usr/bin/perl -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# changes color (hue) to grayscale

if( $h || $help ) {
  print "Usage: hh2g [ -h | -help ] [ -i | -invert ] [ -color ] file\n";
  exit;
}

if( $i || $invert ) {
  $invert = 1;
}

while(<>) {

  if( $color && /^(.+) G$/ ) {
	if( $1 == 0. ) {		#do not transform black
	  print;
	} else {
	  $hsb = &invert( $1 );
	  print "$hsb HH\n";
	}
  } elsif( /^(.+) HH$/ ) {
	$gray = &invert( $1 );
	print "$gray G\n";
  } else {
	print;
  }
}

sub invert {

  my $c = $_[0];

  if( $invert ) {
    $c = 1. - $c;
  }

  return $c;
}
  
