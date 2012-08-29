#!/usr/bin/perl -s
#
# adapted from
# http://www.math.harvard.edu/computing/ps/index.html
#
#------------------------------------------------------------

$text = "Draft" unless defined $text;		# text to write
$fontsize = 72 unless defined $fontsize;	# fontsize of text
$outline = 0 unless defined $outline;		# write text as outline
$shade = 0 unless defined $shade;		# write text shaded
$above = 0 unless defined $above;		# write text above plot

unless( defined $gray ) {
  if( $shade ) {
    $gray = 0.80;
  } else {
    $gray = 0.95;
  }
}

if( $h or $help ) {
  FullUsage();
  exit 0;
} elsif( scalar @ARGV == 0 ) {
  Usage();
  exit 0;
}

#---------------------------------------------------------
# do not change anything below here
#---------------------------------------------------------

$dist = $fontsize + 8;

$start = 1.;
$end = $gray;
$inc = ($end - $start) / 20;

$flag = 0;
while (<>) {
    if ( not $above and /^%%Page:/ ) {
	restore();
	write_watermark();
    } elsif ( $above and /^%PsEndPage/ ) {
	restore();
	write_watermark();
    } else {
        print;
    }
}

sub write_watermark {

        print $_;
        print "% watermark start\n";
        print "gsave\n";
        print "$gray setgray\n";
	print "/printText { 0 0 moveto show } def\n";
        print "/Helvetica-Bold findfont $fontsize scalefont setfont\n";
        print "$dist $dist 800 {\n";
        print "306 exch moveto\n";
        print "($text)\n";
        print "dup stringwidth pop 2 div neg 0 rmoveto\n";
	if( $outline ) {
          print "false charpath $outline setlinewidth stroke\n";
	} elsif( $shade ) {
	  print "gsave\n";
	  print "currentpoint translate\n";
          print "  $start $inc $end {\n";
          print "  setgray -1 0.5 translate\n";
          print "  dup printText\n";
          print "  } for\n";
          print "1 setgray printText\n";
          print "grestore\n";
	} else {
          print "show\n";
	}
        print "} for\n";
        print "grestore\n";
        print "% watermark end\n";
}

sub restore {

        if ($flag) {
            print "grestore\n";
            print "% ggu grestore from flag = 1\n";	#can be deleted
        }
        $flag = 1;
}

sub Usage {
  print STDERR "Usage: watermark.pl [-h|-help] [-options] ps-file\n";
}

sub FullUsage {
  Usage();
  print STDERR "  Inserts watermark text into PostScript file\n";
  print STDERR "  options:\n";
  print STDERR "    -text=string  insert string as watermark (default Draft)\n";
  print STDERR "    -fontsize=#   use # as fontsize (default 72)\n";
  print STDERR "    -above        put text above plot\n";
  print STDERR "    -outline=#    outline text with linewidth # (default 0)\n";
  print STDERR "    -shade        shade text\n";
  print STDERR "    -gray=#       use gray color # (default 0.95,";
  print STDERR " for shade 0.80)\n";
}

