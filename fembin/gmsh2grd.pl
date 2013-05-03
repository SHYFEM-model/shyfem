#!/usr/bin/perl
#
# Converts a mesh in GMSH format (MEDIT program by Inria) into GRD format.
#
# In GMSH use the option to save in *.mesh format and save only
# the physical groups (you must have only 1 surface as physical group!)
#
#-----------------------------------------------------------------------

use Math::Trig;

#ReadOld();
ReadNew();

#---------------------------------------------------

sub ReadNew {		# use this for version 2.2 and higher

  my %hash = ();

  while(<>) {

    if( /^\$Nodes/ ) {
	$totnodes=<>; chomp($totnodes);
	for ($k=0; $k < $totnodes; $k++) {
	    $line=<>; chomp($line);
	    my ($n,$x,$y,$z) =split(/\s+/,$line);
	    ($long,$latg)=&ConvertCoord($x,$y,$z); # (x,y,z)=>(lon,lat)
	    print "1 $n 0 $long $latg\n";
	    $hash{$n} = 1;
	}
    } elsif( /^\$Elements/ ) {
	$totel=<>; chomp($totel);
	for ($ie=0; $ie < $totel; $ie++) {
	    $line=<>; chomp($line);
	    @lline=split(/\s+/,$line);
	    $iel = $lline[0];
	    $ity = $lline[1];
	    if( $ity == 2 ) {
	      print "2 $iel 0 3 $lline[5] $lline[6] $lline[7]\n";
	      $hash{$lline[5]}++;
	      $hash{$lline[6]}++;
	      $hash{$lline[7]}++;
	    }
	}
    }
  }

  foreach (keys %hash) {
    my $c = $hash{$_};
    if( $c <= 1 ) {
      print STDERR "Warning: node $_ has not been used in elements\n";
    }
  }
}

#---------------------------------------------------

sub ReadOld {

  while (<>) {
	if ($_=~/Vertices/) {			# Reads nodes
		$totnodes=<>; chomp($totnodes);
		for ($k=0; $k < $totnodes; $k++) {
		    $line=<>; chomp($line);
		    @lline=split(' ',$line);
		    $kk=$k+1;

		    $x=$lline[0]; $y=$lline[1]; $z=$lline[2];
		    ($long,$latg)=&ConvertCoord($x,$y,$z); # (x,y,z)=>(lon,lat)

		    print "1 $kk 0 $long $latg\n";
		}
	} elsif ($_=~/Triangles/) {		# Reads elements
		$totel=<>; chomp($totel);
		for ($ie=0; $ie < $totel; $ie++) {
		    $line=<>; chomp($line);
		    @lline=split(' ',$line);
		    $iie=$ie+1;

		    #Change if it is not in anticlockwise sense:
		    #print "2 $iie 0 3 $lline[0] $lline[2] $lline[1]\n";
		    print "2 $iie 0 3 $lline[0] $lline[1] $lline[2]\n";
		}
	}
  }
}

#---------------------------------------------------

sub ConvertCoord {

  # Tranform from Cartesian coordinates:

  my ($x,$y,$z) = @_;
  #$x=$_[0]; $y=$_[1]; $z=$_[2];

  my $lon=atan($y/$x);
  my $lat=asin( $z/(sqrt($x**2+$y**2+$z**2)) );

  my $long=($lon / 3.14159265358979)*180;
  my $latg=($lat / 3.14159265358979)*180;

  return($long,$latg);

}

#---------------------------------------------------

