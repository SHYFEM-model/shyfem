#!/usr/bin/perl -s
#
# converts kml (google earth) data to grd format
#
# command line options: -close     close lines
#
#------------------------------------------------------------

$in_coords = 0;
$nnode = 0;
$nline = 0;

while(<>) {

  if( /\<coordinates\>/ ) {
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

