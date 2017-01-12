#!/usr/bin/perl
#
# Converts a mesh in *.msh format # in *.grd format.
#

#define libraries
use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");
use grd;
use Math::Trig;

#define parameters
$grid = new grd;
$type = -1;
$oldtype = -1;

#loop over file lines

while (<>){
 if ($_=~/\$Nodes/){			# Reads nodes
  $totnodes=<>; chomp($totnodes);
  for ($k=0; $k < $totnodes; $k++) {
 	$line=<>; chomp($line);
 	@lline=split(' ',$line);
 	$kk=$k+1;
 	$nn=$lline[0]; $x=$lline[1]; $y=$lline[2]; $z=$lline[3];
 	#($long,$latg)=&ConvertCoord($x,$y,$z);	# Convert coordinates (x,y,z)=>(lon,lat)
 	#print "1 $kk 0 $long $latg\n";
        @ntem = ();
        $ntem[1] = $nn;
        $ntem[2] = 0;
        $ntem[3] = $x;
        $ntem[4] = $y;
   	$grid->insert_node(\@ntem);
   }
 } elsif ($_=~/\$Elements/){		# Reads elements
   $totel=<>; chomp($totel);
   $iie = 0;
   for ($ie=0; $ie < $totel; $ie++) {
 	$line=<>; chomp($line);
 	@lline=split(' ',$line);
 	#Change if it is not in anticlockwise sense:
 	#print "2 $iie 0 3 $lline[0] $lline[2] $lline[1]\n";
 	#print "2 $iie 0 3 $lline[0] $lline[1] $lline[2]\n";
 	#print "2 $iie $lline[3] 3 $lline[0] $lline[1] $lline[2]\n";
        if ($lline[1] eq 2) {
 	  $type++ if ($lline[4] ne $oldtype);
 	  $oldtype = $lline[4];
 	  $iie=$iie+1;
          @ntem = ();
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

$grid->delete_unused();			#delete unused node
#$grid->set_preserve_order($::preserve);
$grid->writegrd("gsmh_mesh.grd");	#write grd mesh

#################################
#Routine not used
sub ConvertCoord{
$x=$_[0]; $y=$_[1]; $z=$_[2];

# Tranform from Cartesian coordinates:
$lon=atan($y/$x);
$lat=asin( $z/(sqrt($x**2+$y**2+$z**2)) );

$long=($lon / 3.14159265358979)*180;
$latg=($lat / 3.14159265358979)*180;
return($long,$latg);
}
