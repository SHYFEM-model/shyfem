#!/usr/bin/perl -s
#
# inserts eps file into another ps/eps file
#
# (C) Georg Umgiesser 2001
#
# Usage: gpsinsert.pl [-relative] "x1 y1 x2 y2" eps-insert-file ps-file
#
#	-relative		bounding box is relative [0-1]
#	"x1 y1 x2 y2"		area where to insert (in point units)
#	eps-insert-file		eps file to be inserted
#	ps-file			ps-file into which eps-insert-file is inserted
#
# result is written to stdout

$where = shift;
$insertfile = shift;
$psfile = shift;

open(EPS,"<$insertfile");
@insertfile = <EPS>;
close(EPS);

open(PS,"<$psfile");
@psfile = <PS>;
close(PS);

$bbps = setbox(\@psfile);
$bb = setbox(\@insertfile);
print STDERR "Bounding Box PS: $bbps\n";
print STDERR "Bounding Box: $bb - Landscape: $landscape\n";

($wx1,$wy1,$wx2,$wy2) = &getbox($where);
($bx1,$by1,$bx2,$by2) = &getbox($bb);
($px1,$py1,$px2,$py2) = &getbox($bbps);

if( $landscape ) {
  ($wx1,$wy1) = &swap($wx1,$wy1);
  ($wx2,$wy2) = &swap($wx2,$wy2);
}

if( $relative ) {
  print STDERR "relative coordinates ... converting\n";
  $wx1 = $px1 + $wx1 * ($px2-$px1);
  $wx2 = $px1 + $wx2 * ($px2-$px1);
  $wy1 = $py1 + $wy1 * ($py2-$py1);
  $wy2 = $py1 + $wy2 * ($py2-$py1);
}

$dwx = $wx2 - $wx1;
$dwy = $wy2 - $wy1;
$dbx = $bx2 - $bx1;
$dby = $by2 - $by1;
$sx = $dwx / $dbx;
$sy = $dwy / $dby;

# rectify scale if where-box has negative width

if( $dwx <= 0 && $dwy <= 0 ) { 
  die "Impossible to insert: $where\n";
} elsif( $dwx <= 0 ) {
  $sx = $sy;
  $dwx = $sx * $dbx;
  $wx2 = $wx1 + $dwx;
} elsif( $dwy <= 0 ) {
  $sy = $sx;
  $dwy = $sy * $dby;
  $wy2 = $wy1 + $dwy;
}

if( $landscape ) {
  ($wx1,$wy1) = &swap($wx1,$wy1);
  ($wx2,$wy2) = &swap($wx2,$wy2);
  ($dwx,$dwy) = &swap($dwx,$dwy);
}

print STDERR "$wx1,$wy1,$wx2,$wy2,$dwx,$dwy\n";

foreach (@psfile) {

  if( /^%%Creator: (\S+)/ ) {
    $creator = $1;
    if( $creator eq "psgraph" ) { 
	$ispsgraph = 1;
    } elsif( $creator eq "gnuplot" ) { 
	die "Cannot insert into gnuplot eps\n";
    } else {
	die "Unknown creator $creator - cannot insert\n";
    }
    print STDERR "Creator: $creator    $ispsgraph\n";
  }

  if( $ispsgraph ) {
    if( /showpage$/ ) {
      print "% start inserting eps file (ggu)\n";
      print "BeginEPSF\n";
      print "$wx1 $wy1 translate\n";
      print "-90 rotate\n" if $landscape;
      print "$sx $sy scale\n";
      print "-$dbx 0 translate\n" if $landscape;
      print "-$bx1 -$by1 translate\n";
      print "1 setgray\n";	#1, 0 test
      print "$bx1 $by1 $dbx $dby Rect fill\n";
      print "0 setgray\n";
      print "$bx1 $by1 $dbx $dby Rect clip newpath\n";
      print "\n";
      print "%%BeginDocument: $insertfile\n";
      foreach $line (@insertfile) {
	  print $line;
      }
      print "%%EndDocument\n";
      print "\n";
      print "EndEPSF\n";
    }
  }

  print;
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
  $landscape = 0;

  foreach my $line (@$file) {
    if( $line =~ /%%BoundingBox:\s+(.*)/ ) {
      if( $1 eq "(atend)" ) {
	;
      } else {
	$bb = $1;
      }
    }
    if( $line =~ /^%%Orientation: Landscape/ ) {
      $landscape = 1;
    }
  }

  die "No bounding box found in insert file\n" unless $bb;

  return $bb;
}

###################################################################

sub swap {
  my ($x,$y) = @_;
  return ($y,$x);
}

###################################################################

