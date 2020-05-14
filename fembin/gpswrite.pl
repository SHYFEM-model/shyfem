#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# writes extra text to ps/eps file
#
# uses absolute coordinates

$text = "(georg) 300 300 50 AW";
$text = shift;

while(<>) {

  if( /^%%Creator: gnuplot/ ) { $isgnuplot = 1; }
  if( /^%%Creator: psgraph/ ) { $ispsgraph = 1; }

  if( $isgnuplot ) {
    if( /^\/HVSF.*def$/ ) { $hasprolog = 1; }
    if( /^end$/ ) {
      $end++;
      unless( $hasprolog ) {
        print "/HVSF { (Helvetica) findfont exch scalefont setfont } def\n";
        print "/AM { 2 copy 2 copy translate 90 rotate\n";
        print "  neg exch neg exch translate moveto } def % Absolute Move\n";
        print "/AW { gsave HVSF AM show grestore } def % (text) x y points\n";
        $hasprolog = 1;
      }
      if( $end == 2 ) {
        print "% extra test inserted (ggu)\n";
        print "$text";
        print "\n";
      }
    }
  }

  if( $ispsgraph ) {
    if( /^\/TRSF.*def$/ ) { $hasprolog = 1; }
    if( /^%%EndResource$/ ) { 	#for old format
      unless( $hasprolog ) {
	print "/TRSF { /Times-Roman FF exch CF SF } bind def\n";
	print "/AW { TRSF M show } bind def % (text) x y points\n";
        $hasprolog = 1;
      }
    }
    if( /^showpage$/ ) {
      print "% extra test inserted (ggu)\n";
      print "$text";
      print "\n";
    }
    if( /^GR showpage$/ ) {	#for old format
      print "GR\n";
      print "% extra test inserted (ggu)\n";
      print "$text";
      print "\n";
      print "showpage\n";
      next;
    }
  }

  print;
}
