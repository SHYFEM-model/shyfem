#!/usr/bin/perl
#
# reads dxf file and extracts shapes and writes to grd file
#
# can read the following shapes:
#
#	AcDbPolyline
#
#----------------------------------------------------------

%blocks = (
		 "SECTION"	=>	"ENDSEC"
		,"TABLE"	=>	"ENDTAB"
		,"BLOCK"	=>	"ENDBLK"
	  );

# 100 AcDbPolyline
# 100 AcDbXrecord

$space = "";
$end = "";

$::nodes = 0;
$::lines = 0;

open(GRD,">dxf.grd") || die "Cannot open dfx.grd file\n";

while(1) {

  my ($key,$value) = get_entity();

  last if $key < 0;

  if( $key == 0 ) {
    if( $value eq $end ) {
      $end = pop(@end);
      $space = pop(@space);
      print STDERR "$space $value\n";
    } elsif( $blocks{$value} ) {
      push(@end,$end);
      $end = $blocks{$value};
      print STDERR "$space $value\n";
      push(@space,$space);
      $space .= "  ";
    } else {
      print STDERR "$space $value\n";
    }
  } elsif( $key == 2 ) {
    print STDERR "$space $value ($key)\n";
  } elsif( $key == 100 ) {
    print STDERR "$space $value ($key)\n";
    if( $value eq "AcDbPolyline" ) {
      print STDERR "$space ********* Polyline.........\n";
      parse_polyline();
    }
  } else {
    $zero = 0;
  }
}

close(GRD);

#----------------------------------------------------------

sub get_entity {

  $_ = <>;
  return -1 unless $_;

  chomp;
  s/\r$//;

  my $key = $_;

  $_ = <>;
  chomp;
  s/\r$//;

  my $value = $_;

  return ($key,$value);
}

sub get_known_entity {

  my $expect = shift;

  my ($key,$value) = get_entity();

  if( $key != $expect ) {
    die "*** cannot parse polyline for $expect: $key\n";
  }

  return $value;
}

#----------------------------------------------------------

sub parse_polyline {

  my @nodes = ();

  my $n = get_known_entity(90);
  my $aux = get_known_entity(70);
  my $aux = get_known_entity(43);
  my $depth = get_known_entity(38);

  print GRD "\n";

  for( my $i=0; $i<$n; $i++ ) {
    my $x = get_known_entity(10);
    my $y = get_known_entity(20);
    $::nodes++;
    push(@nodes,$::nodes);
    print GRD "1 $::nodes 0 $x $y $depth\n";
  }

  $n++;			#close line
  push(@nodes,$nodes[0]);

  print GRD "\n";

  $::lines++;
  print GRD "3 $::lines 0 $n\n";

  my $i = 0;
  foreach my $node (@nodes) {
    $i++;
    print GRD " $node";
    print GRD "\n" if $i%10 == 0;
  }
  print GRD "\n";
  
}








