#!/usr/bin/perl -w
# Convert boundary line in the format for GMSH
# BSpline is the default line type, anyway it could be changed manually 
# in the file. Closed lines need to be defined as BSpline.

use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");

use grd;
use strict;

my $grid = new grd;
my $file = $ARGV[0];

$grid->readgrd($file);

my $lines = $grid->get_lines();
my $nodes = $grid->get_nodes();

print "// --------------------------------------------------//\n";
print "// NODE LIST //\n";
print "// --------------------------------------------------//\n";

# Print node list #
foreach my $node (sort {$a <=> $b} values %$nodes) {
    my $x = $node->{x};
    my $y = $node->{y};
    my $nnum = $node->{number};
    print "Point($nnum) = {$x, $y, 0, 0.0};\n";
}
print "\n";

print "// --------------------------------------------------//\n";
print "// LINE LIST //\n";
print "// --------------------------------------------------//\n";
print "\n";

# Print line list #
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
  #print "Line($lnum) = {";
  print "BSpline($lnum) = {";
  print join( ',', @nlist );
  print "};\n";
}

print "\n";
print "// --------------------------------------------------//\n";
print "// LINE LOOPS AND PLANE SURFACES //\n";
print "// --------------------------------------------------//\n";


print "\n";
print "// --------------------------------------------------//\n";
print "// SIZE FIELDS //\n";
print "// --------------------------------------------------//\n";

print "\n";
print "// --------------------------------------------------//\n";
print "// MESHING PROPERTIES //\n";
print "// --------------------------------------------------//\n";
print "Geometry.LineNumbers = 1;\n";
print "Mesh.Algorithm = 6;\n";
print "Mesh.ChacoHypercubeDim = 0;\n";
print "Mesh.ChacoMeshDim1 = 1;\n";
print "Mesh.CharacteristicLengthExtendFromBoundary = 0;\n";
print "Mesh.CharacteristicLengthFromPoints = 0;\n";
print "Mesh.RemeshParametrization = 0;\n";
print "Mesh.SurfaceFaces = 1;\n";
print "\n";
