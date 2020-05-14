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
# converts kml (google earth) data to grd format
#
# command line options: -close           close lines
# command line options: -points_only     only write points
#
#------------------------------------------------------------

$in_coords = 0;
$nnode = 0;
$nline = 0;

fullusage() if $h or $help;
usage() unless $ARGV[0];

while(<>) {

  if ( /\<coordinates\>(.+?)\<\/coordinates\>/gi ) {
    $coords = $1;
    parse_coords($coords);
  } elsif( /\<coordinates\>/ ) {
    $coords = <>;

    $line = <>;
    unless( $line =~ /\<\/coordinates\>/ ) {
	die "expecting closure of coordinates... $ARGV\n";
    }

    parse_coords($coords);
  }
}

print STDERR "total nodes: $nnode  total lines: $nline\n";

#------------------------------------------------------

sub parse_coords {

  my ($line) = @_;

  $line =~ s/^\s+//;

  my @f = split(/\s+/,$line);
  my @nodes = ();
  my ($x,$y,$z);
  my ($x0,$y0,$z0);
  my $nelim = 0;

  my $n = @f;
  my $first = 1;

  print "\n";
  foreach my $coord (@f) {
    ($x,$y,$z) = split(/,/,$coord);
    if( $first ) {
      ($x0,$y0,$z0) = ($x,$y,$z);
      $first = 0;
    } elsif( $x == $x0 and $y == $y0 and $z == $z0 ) {	# eliminate
      $nelim = $nnode;
      next;
    }
    $nnode++;
    push(@nodes,$nnode);
    #print STDERR "$nnode   $x $y $z\n";
    print "1 $nnode 0 $x $y $z\n";
  }

  return if $points_only;

  if( $close or $nelim == $nnode ) {
    #print STDERR "closing line...\n";
    push(@nodes,$nodes[0]);
  }

  $nline++;
  $n = @nodes;
  print STDERR "file $ARGV  line $nline read with nodes $n\n";

  print "\n";
  print "3 $nline 0 $n\n";

  my $i = 0;
  foreach my $node (@nodes) {
    $i++;
    print " $node";
    print "\n" if $i%10 == 0;
  }

  print "\n";
}

#------------------------------------------------------

sub fullusage {
  print "Usage: kml2grd.pl [-h|-help] [-options] file(s)\n";
  print "  Transforms KML files (Google Earth) to GRD format\n";
  print "    -h|-help         this help message\n";
  print "    -close           close lines\n";
  print "    -points_only     only write points, make no lines\n";
  print "    file(s)          file(s) in KML (Google Earth) format\n";
  exit 0;
}
sub usage {
  print "Usage: kml2grd.pl [-h|-help] [-options] file(s)\n";
  print "  Transforms KML files (Google Earth) to GRD format\n";
  exit 0;
}

#------------------------------------------------------

