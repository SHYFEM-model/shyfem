#!/usr/bin/perl -ws
#
# Converts a GMSH mesh in *.msh format in *.grd format.
# Check the created .grd mesh with the command grid -k gsmh_msh.grd
# for clockwise elements and node connections 
#
#--------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
use Math::Trig;
use strict;

$::help = 1 if $::help;
$::h = 1 if $::h;
$::geo = 0 unless $::geo;
$::line = 0 unless $::line;

if( $::h or $::help ) {
  Usage();
    exit 0;
} elsif( not $ARGV[0] ) {
  Usage();
    exit 1;
  print STDERR "No msh file given\n";
}

#-------------------------------------------------

#define parameters
my $grid = new grd;
my $type = 0;
my $oldtype = 0;

#-------------------------------------------------
# Loop over file lines
#-------------------------------------------------

while (<>){
 if ($_=~/\$Nodes/){			# Reads nodes
  my $totnodes=<>; chomp($totnodes);
  for (my $k=0; $k < $totnodes; $k++) {
 	my $line=<>; chomp($line);
 	my @lline=split(' ',$line);
 	my $nn=$lline[0]; my $x=$lline[1]; my $y=$lline[2]; my $z=$lline[3];
        my @ntem = ();
        $ntem[1] = $nn;
        $ntem[2] = 0;
        if( $::geo ) {
  	  my ($long,$latg)=&ConvertCoord($x,$y,$z);
          $ntem[3] = $long;
          $ntem[4] = $latg;
        } else {
          $ntem[3] = $x;
          $ntem[4] = $y;
        }
   	$grid->insert_node(\@ntem);
   }
 } elsif ($_=~/\$Elements/){		# Reads elements
   my $totel=<>; chomp($totel);
   my $iie = 0;
   for (my $ie=0; $ie < $totel; $ie++) {
 	my $line=<>; chomp($line);
 	my @lline=split(' ',$line);
 	#Change if it is not in anticlockwise sense:
 	#print "2 $iie 0 3 $lline[0] $lline[2] $lline[1]\n";
 	#print "2 $iie 0 3 $lline[0] $lline[1] $lline[2]\n";
 	#print "2 $iie $lline[3] 3 $lline[0] $lline[1] $lline[2]\n";
        if ($lline[1] eq 2) {
 	  $type++ if ($lline[4] ne $oldtype);
 	  $oldtype = $lline[4];
 	  $iie=$iie+1;
          my @ntem = ();
          $ntem[1] = $iie;
          $ntem[2] = $type;
          $ntem[3] = 3;
          $ntem[4] = $lline[5];
          $ntem[5] = $lline[6];
          $ntem[6] = $lline[7];
          $grid->insert_elem(\@ntem);
        }
   }
 }
}
#-------------------------------------------------
# Write grid to fle gsmh_msh.grd
#-------------------------------------------------
$grid->delete_unused();			#delete unused node
$grid->writegrd("gsmh_msh.grd");

#--------------------------------------
#--------------------------------------
# Tranform from stereographic coordinates:
sub ConvertCoord{
  my $x=$_[0]; my $y=$_[1]; my $z=$_[2];
  my $lon=atan($y/$x);
  my $lat=asin( $z/(sqrt($x**2+$y**2+$z**2)) );
  my $deg2rad = 3.14159265358979 / 180.;
  my $long=$lon / $deg2rad;
  my $latg=$lat / $deg2rad;
  return($long,$latg);
}

#--------------------------------------

sub Usage {
  print STDERR "                                    \n";
  print STDERR "Usage: gmsh2grd.pl [-options] msh_file\n";
  print STDERR "  -h|-help  this help screen        \n";
  print STDERR "  -geo:     coordinates in lat/lon  \n";
  print STDERR "                                    \n";
  print STDERR "Additional info:                    \n";
  print STDERR " This script creates the file gsmh_msh.grd.   \n";
  print STDERR " Check the mesh with the command grid -k gsmh_msh.grd\n";
  print STDERR " for clockwise elements and node connections. \n";
}
