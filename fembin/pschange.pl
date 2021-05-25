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
# changes color, line thickness, text size in ps files created by gnuplot
#
#------------------------------------------------------------

$size_old = 140;		# original text size in ps
$size_new_general = 160;	# desired general text size (axis, title, etc..)
$size_new_plot = 180;		# desired text size of plot labels
$width = 0;			# line width (original is 1)
$change_color = 0;		# change line from black to color

# changes are done only if value > 0

#------------------------------------------------------------

$plot = 0;
$debug = 0;

while(<>) {

  if( /^% Begin plot #(\w+)/ ) {
    $plot = $1;
    print STDERR "starting reading plot $plot\n";
  }

  if( /^% End plot #(\w+)/ ) {
    die "error in numbering of plots: $1 $plot\n" if( $1 != $plot );
    print STDERR "finished reading plot $plot\n";
    $plot = 0;
  }

  if( $plot ) {					# ------- inside plot area
    if( /^LCb/ and $change_color ) {		# change line color
      $color = "LC" . $plot;			# use color $plot
      print STDERR "  $_" if $debug;
      s/LCb/$color/;
      print STDERR "  $_" if $debug;
    }
    if( / UL/ and $width ) {			# change line width
      print STDERR "  $_" if $debug;
      $_ = "$width UL\n";
      print STDERR "  $_" if $debug;
    }
    if( /Helvetica/ and $size_new_plot ) {	# change text size in label
      print STDERR "  $_" if $debug;
      s/$size_old/$size_new_plot/;
      print STDERR "  $_" if $debug;
    }
  } else {					# ------- outside plot area
    if( /Helvetica/ and $size_new_general ) {	# change general text size
      print STDERR "  $_" if $debug;
      s/$size_old/$size_new_general/;
      print STDERR "  $_" if $debug;
    }
  }

  print;
}

