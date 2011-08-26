#!/usr/bin/perl -s
#
# inserts eps file into another ps/eps file
#
# (C) Georg Umgiesser 2010-2011
#
# Usage: gpsinsert.pl [options] "x1 y1 x2 y2" eps-insert-file ps-file
#
#	options:
#	  -relative		bounding box is relative [0-1]
#	  -reggrid		inserts reggrid into file for orientation
#	"x1 y1 x2 y2"		bounding box where to insert (in point units)
#	eps-insert-file		eps file to be inserted
#	ps-file			ps-file into which eps-insert-file is inserted
#
# result is written to stdout
# if one of x2/y2 is zero the imported file conserves its scale factor
#
# version 2.1
#
# 07.10.2010	ggu	allow for more than one insert
# 25.08.2011	ggu	insert in all ps files, include code directly
#
# still to test:
#   - insert in any ps file
#   - insert landscape eps file
#
#-------------------------------------------------------------------

$where = shift;
$insertfile = shift;
$psfile = shift;

open(EPS,"<$insertfile");
@insertfile = <EPS>;
close(EPS);

open(PS,"<$psfile");
@psfile = <PS>;
close(PS);

($bbps,$landscape_ps) = setbox(\@psfile);
($bbeps,$landscape_eps) = setbox(\@insertfile);

print STDERR "Bounding Box PS: $bbps - Landscape: $landscape_ps\n";
print STDERR "Bounding Box EPS: $bbeps - Landscape: $landscape_eps\n";

($wx1,$wy1,$wx2,$wy2) = getbox($where);
($bx1,$by1,$bx2,$by2) = getbox($bbeps);
($px1,$py1,$px2,$py2) = getbox($bbps);

#----------------------------------------------- handle landscape eps file

if( $landscape_eps ) {			# still to be checked
  ($bx1,$by1) = swap($bx1,$by1);
  ($bx2,$by2) = swap($bx2,$by2);
}

#----------------------------------------------- convert relative to absolute

if( $landscape_ps ) {
  ($px1,$py1) = swap($px1,$py1);	# invert just for rectyfying
  ($px2,$py2) = swap($px2,$py2);
}

if( $relative ) {
  print STDERR "relative coordinates ... converting\n";
  $wx1 = $px1 + $wx1 * ($px2 - $px1);
  $wx2 = $px1 + $wx2 * ($px2 - $px1);
  $wy1 = $py1 + $wy1 * ($py2 - $py1);
  $wy2 = $py1 + $wy2 * ($py2 - $py1);
}

$dbx = $bx2 - $bx1;
$dby = $by2 - $by1;
$dwx = $wx2 - $wx1;
$dwy = $wy2 - $wy1;
$sx = $dwx / $dbx;
$sy = $dwy / $dby;

print STDERR "before rectifying: $wx1 $wy1   $wx2 $wy2   $dwx $dwy\n";

#----------------------------------------------- rectify scale for dw<=0

if( $dwx <= 0 and $dwy <= 0 ) { 
  die "Impossible to insert: $where ... at least one dimension must be > 0\n";
} elsif( $dwx <= 0 ) {
  $sx = $sy;
  $dwx = $sx * $dbx;
  $wx2 = $wx1 + $dwx;
} elsif( $dwy <= 0 ) {
  $sy = $sx;
  $dwy = $sy * $dby;
  $wy2 = $wy1 + $dwy;
}

print STDERR "insert window: $wx1 $wy1   $wx2 $wy2   $dwx $dwy\n";
print STDERR "insert box:    $bx1 $by1   $bx2 $by2   $dbx $dby\n";

if( $landscape_ps ) {			# transform everything to physical
  ($px1,$py1) = swap($px1,$py1);
  ($px2,$py2) = swap($px2,$py2);
  ($wx1,$wy1) = land2port($wx1,$wy1);
  ($wx2,$wy2) = land2port($wx2,$wy2);
  ($dwx,$dwy) = swap($dwx,$dwy);
  print STDERR "insert window: $wx1 $wy1   $wx2 $wy2   $dwx $dwy\n";
  print STDERR "insert box:    $bx1 $by1   $bx2 $by2   $dbx $dby\n";
}

#----------------------------------------------- start inserting

$creator = "";
$in_epsf = 0;				# only insert in original ps file
$code_inserted = 0;
$page_inserted = 0;

foreach (@psfile) {

  if( not $creator and /^%%Creator: (\S+)/ ) {
    $creator = $1;
    if( $creator eq "psgraph" ) { 
	;
    } elsif( $creator eq "gnuplot" ) { 
	;
    } else {
	print STDERR "Unknown creator $creator - trying to insert\n";
    }
    print STDERR "Creator: $creator\n";
  }

  if( /^%%EndProlog/ and not $in_epsf ) {	# insert code
    insert_EPSF();
    $code_inserted = 1;
  }

  $in_epsf = 1 if /ggu start-eps/;	# inside EPSF section
  $in_epsf = 0 if /ggu end-eps/;	# outside EPSF section

  if( not $in_epsf ) {
    if( /showpage$/ and $code_inserted ) {	# only if code was inserted
      $page_inserted = 1;
      print "% start inserting eps file (ggu start-eps)\n";

      print "$px1 $py1 $px2 $py2 RegGrid\n" if $reggrid;	# only for test

      print "BeginEPSF\n";

      print "$wx1 $wy1 translate\n";
      print "90 rotate\n" if $landscape_ps;
      print "$sx $sy scale\n";

      # next two lines not sure
      #print "$dbx 0 translate\n" if $landscape_ps;
      print "-$bx1 -$by1 translate\n";

      print "1 setgray\n";			# 1=white, 0=black (for test)
      print "$bx1 $by1 $dbx $dby Rect fill\n";
      print "0 setgray\n";
      print "$bx1 $by1 $dbx $dby Rect clip newpath\n";
      print "\n";

      print "% start of external eps file (ggu)\n";
      print "%%BeginDocument: $insertfile\n";
      foreach $line (@insertfile) {
	  print $line;
      }
      print "%%EndDocument\n";
      print "% end of external eps file (ggu end-eps)\n";
      print "\n";

      print "EndEPSF\n";
    }
  }

  print;
}

#----------------------------------------------- error messages

if( not $code_inserted or not $page_inserted ) {
  if( not $code_inserted ) {
    print STDERR "Could not insert code... unknown format of ps file\n";
  } elsif( not $page_inserted ) {
    print STDERR "Could not insert eps-file... no showpage in ps file?\n";
  }
  print STDERR "Creator: $creator\n";
}

###################################################################

sub land2port {

  my ($x,$y) = @_;

  my $aux = $y - $px1;
  $y = $x;
  $x = $px2 - $aux;

  return ($x,$y);
}

###################################################################

sub getbox {

  $_ = shift;

  s/^\s+//;

  my @f = split;
  my $n = @f;

  die "box must have 4 values (x1,y1,x2,y2): $_\n" if $n != 4;

  return @f;
}

###################################################################

sub setbox {

  my $file = shift;

  my $line;
  my $bb = "";
  my $in_epsf = 0;
  my $landscape = 0;

  foreach my $line (@$file) {
    $in_epsf = 1 if $line =~ /ggu start-eps/;	# inside EPSF section
    $in_epsf = 0 if $line =~ /ggu end-eps/;	# outside EPSF section

    if( not $in_epsf and $line =~ /%%BoundingBox:\s+(.*)/ ) {
      if( $1 eq "(atend)" ) {
	;
      } else {
	$bb = $1;
      }
    }
    if( not $in_epsf and $line =~ /^%%Orientation: Landscape/ ) {
      $landscape = 1;
    }
  }

  die "No bounding box found in insert file\n" unless $bb;

  return ($bb,$landscape);
}

###################################################################

sub swap {
  my ($x,$y) = @_;
  return ($y,$x);
}

###################################################################

sub insert_EPSF {	# inserts code to include eps files
print <<EOI
%%BeginResource: eps_procset_new
/BeginEPSF {
  /b4_Inc_state save def
  /dict_count countdictstack def
  /op_count count 1 sub def
  userdict begin
  /showpage { } def
  0 setgray 0 setlinecap
  1 setlinewidth 0 setlinejoin
  10 setmiterlimit [ ] 0 setdash newpath
  languagelevel where
  { pop languagelevel
  1 ne
    {false setstrokeadjust false setoverprint
    } if
  } if
} bind def
/EndEPSF {
  count op_count sub {pop} repeat
  countdictstack dict_count sub {end} repeat
  b4_Inc_state restore
} bind def
/Rect { % llx lly w h
  4 2 roll moveto 1 index 0 rlineto
  0 exch rlineto neg 0 rlineto closepath
} bind def
/ClipRect { % llx lly w h
  Rect clip newpath
} bind def
/RegGrid { % bx1 by1 bx2 by2
/by2 exch def /bx2 exch def /by1 exch def /bx1 exch def
/dbx { bx2 bx1 sub 10 div } def
/dby { by2 by1 sub 10 div } def
0.5 setgray
newpath bx1 11 { dup by1 moveto dup by2 lineto dbx add } repeat stroke
newpath by1 11 { dup bx1 exch moveto dup bx2 exch lineto dby add } repeat stroke
0 setgray
} bind def
%%EndResource
EOI
;
}
