#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# Convert boundary line in the format for GMSH
# BSpline is the default line type, anyway it could be changed by
# setting the option -Line. Closed lines need to be defined as 
# BSpline.
#
#-----------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use grd;
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
  print STDERR "No line given\n";
}

#-------------------------------------------------

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $lines = $grid->get_lines();
my $nodes = $grid->get_nodes();

my $gmshf = "gmsh_line.geo";
open(GMSH, ">$gmshf");

#-------------------------------------------------
# Print node list 
#-------------------------------------------------

print GMSH "// --------------------------------------------------//\n";
print GMSH "// NODE LIST //\n";
print GMSH "// --------------------------------------------------//\n";

if( $::geo ) {
  my $kmax =  $grid->{nnmax};
  my $ksp = $kmax + 1;
  my $knp = $ksp + 1;
  print GMSH "// Lon/Lat converted to stereographic coordinates.  // \n";
  print GMSH "// Sphere with center at Point($ksp) and           // \n"; 
  print GMSH "// North Pole at Point($knp)                       //\n";
  print GMSH "IS = news;\n";
  print GMSH "Point($ksp) = {0.0,0.0,0.0};         //The center of the sphere\n";
  print GMSH "Point($knp) = {0.0,0.0,6.37101e+06}; //The North Pole\n";
  print GMSH "PolarSphere(1) = {$ksp , $knp};\n"; 
  print GMSH "\n";
}

foreach my $node (sort {$a <=> $b} values %$nodes) {
    my $nnum = $node->{number};
    my $x = $node->{x};
    my $y = $node->{y};
    if( $::geo ) {
      my($xs,$ys)=&ConvertCoord($x,$y);
      print GMSH "Point($nnum) = {$xs, $ys, 0 };\n";
    } else {
      print GMSH "Point($nnum) = {$x, $y, 0, 0.0};\n";
    }
}
print GMSH "\n";

#-------------------------------------------------
# Print line list 
#-------------------------------------------------

print GMSH "// --------------------------------------------------//\n";
print GMSH "// LINE LIST //\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "\n";

foreach my $line (values %$lines) {
  my $nvert = $line->{nvert};
  my $lnum = $line->{number};
  my $vert = $line->{vert};
  my @nlist = ();
  foreach my $nnode (@$vert) {
    my $node = $grid->get_node($nnode);
    my $nnum = $node->{number};
    push(@nlist,$nnum);
  }
  if( $::line ) {
    print GMSH "Line($lnum) = {";
  } else {
    print GMSH "BSpline($lnum) = {";
  }
  print GMSH join( ',', @nlist );
  print GMSH "};\n";
}

#-------------------------------------------------
# Print header of section LINE LOOPS AND SURFACES 
#-------------------------------------------------

print GMSH "\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "// LINE LOOPS AND SURFACES //\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "// The Line Loops and the Plain Surface defined with //\n";
print GMSH "// the GMSH gui must be placed before the Size Field.//\n";
print GMSH "// Remember also to manually add Physical Surface    //\n";
print GMSH "// Example:                                          //\n";
print GMSH "//    Plane Surface(739) = {736, 737, 738};          //\n";
print GMSH "//    Physical Surface(740) = {739};                 //\n";

#-------------------------------------------------
# Print header of section SIZE FIELDS with commented examples 
#-------------------------------------------------

print GMSH "\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "// SIZE FIELDS //\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "// ATTRACTOR to line //\n";
print GMSH "//Field[1] = Attractor;\n";
print GMSH "//Field[1].EdgesList = {2, 7, 8, 6};\n";
print GMSH "//Field[1].NNodesByEdge = 50000;\n";
print GMSH "//Field[2] = Threshold;\n";
print GMSH "//Field[2].DistMax = 3000;\n";
print GMSH "//Field[2].DistMin = 0;\n";
print GMSH "//Field[2].IField = 1;\n";
print GMSH "//Field[2].LcMax = 2000;\n";
print GMSH "//Field[2].LcMin = 2000;\n";
print GMSH "//Field[2].Sigmoid = 1;\n";
print GMSH "//Field[2].StopAtDistMax = 1;\n";
print GMSH "\n";
print GMSH "//MathEval field for defining general resolution over whole system //\n";
print GMSH "//Field[3] = MathEval;\n";
print GMSH "//Field[3].F = \"1300\";\n";
print GMSH "\n";
print GMSH "// MathEval inside specific surface //\n";
print GMSH "//Field[4] = MathEval;\n";
print GMSH "//Field[4].F = \"150\";\n";
print GMSH "//Field[5] = Restrict;\n";
print GMSH "//Field[5].FacesList = {231};\n";
print GMSH "//Field[5].IField = 9;\n";
print GMSH "\n";
print GMSH "// MinAniso to combine all size fields together //\n";
print GMSH "//Field[10] = MinAniso;\n";
print GMSH "//Field[10].FieldsList = {2, 3, 5};\n";
print GMSH "\n";
print GMSH "// Define Background Field //\n";
print GMSH "//Background Field = 10;\n";

#-------------------------------------------------
# Print header of section MESHING PROPERTIES
#-------------------------------------------------

print GMSH "\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "// MESHING PROPERTIES //\n";
print GMSH "// --------------------------------------------------//\n";
print GMSH "Geometry.LineNumbers = 1;\n";
print GMSH "Mesh.Algorithm = 6;\n";
print GMSH "Mesh.ChacoHypercubeDim = 0;\n";
print GMSH "Mesh.ChacoMeshDim1 = 1;\n";
print GMSH "Mesh.CharacteristicLengthExtendFromBoundary = 0;\n";
print GMSH "Mesh.CharacteristicLengthFromPoints = 0;\n";
print GMSH "Mesh.RemeshParametrization = 0;\n";
print GMSH "Mesh.LcIntegrationPrecision = 0.001; // reduce 1d mesh accuracy to speed things up ;\n";
print GMSH "Mesh.SurfaceFaces = 1;\n";
print GMSH "\n";

#-------------------------------------------------
# Convert lon/lat to stereographic coordinates for GMSH
#-------------------------------------------------

sub ConvertCoord{
  my $x=$_[0]; my $y=$_[1];
  my $deg2rad = 3.14159265358979 / 180.;
  my $lon = $x * $deg2rad;
  my $lat = $y * $deg2rad;
  $x = cos($lon)*cos($lat)/(1 + sin($lat));
  $y = cos($lat)*sin($lon)/(1. + sin($lat));
  return($x,$y);
}

#-------------------------------------------------
# Usage
#-------------------------------------------------

sub Usage {
  print STDERR "                                      \n";
  print STDERR "Usage: grd2gmsh.pl [-options] grd_line\n";
  print STDERR "  -h|-help  this help screen          \n";
  print STDERR "  -geo:     coordinates in  lat/lon   \n";
  print STDERR "  -stline:  use straight Line instead of BSpline\n";
  print STDERR "                                      \n";
  print STDERR "Additional info:                      \n";
  print STDERR " This script creates the file gmsh_line.geo.  \n";
  print STDERR " You have to create Plain Surface with GMSH gui. Place them\n";
  print STDERR " before the Size Field and add manually Physical Surface.\n";
}


