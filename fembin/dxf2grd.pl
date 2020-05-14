#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# reads dxf file and extracts shapes and writes to grd file
#
# can read the following shapes:
#
#	AcDbPolyline
#	POLYLINE
#	POINT
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
    if( $value eq "POLYLINE" ) {
      parse_POLYLINE();
    }
    if( $value eq "POINT" ) {
      parse_POINT();
    }
  } elsif( $key == 2 ) {
    print STDERR "$space $value ($key)\n";
  } elsif( $key == 100 ) {
    print STDERR "$space $value ($key)\n";
    if( $value eq "AcDbPolyline" ) {
      parse_AcDbPolyline();
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
    die "*** unexpected key: expected $expect   found $key\n";
  }

  return $value;
}

sub check_known_entity {

  my ($expect,$expval) = @_;

  my $value = get_known_entity($expect);

  if( $value ne $expval ) {
    die "*** unexpected value: expected $exval   found $value\n";
  }
}

#----------------------------------------------------------

sub write_line {

  my ($nodes,$type) = @_;

  $type = 0 unless $type;
  my $n = @$nodes;

  $::lines++;
  print GRD "3 $::lines $type $n\n";

  my $i = 0;
  foreach my $node (@$nodes) {
    $i++;
    print GRD " $node";
    print GRD "\n" if $i%10 == 0;
  }
  print GRD "\n";
}

#----------------------------------------------------------

sub parse_AcDbPolyline {

  #print STDERR "$space ********* AcDbPolyline.........\n";

  my @nodes = ();

  my $n = get_known_entity(90);
  my $aux1 = get_known_entity(70);
  my $aux2 = get_known_entity(43);
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

  write_line(\@nodes);
}

sub parse_POLYLINE {

  # when writing nodes we have to defer writing until we know we
  # are at the last point... this is written only if not equal to the first

  my @nodes = ();
  my ($x,$y,$z);
  my ($x0,$y0);

  my $aux1 = get_known_entity(8);
  my $aux2 = get_known_entity(66);

  print GRD "\n";

  while(1) {
    my $value = get_known_entity(0);
    last if $value eq "SEQEND";
    die "*** expecting different value: $value\n" if $value ne "VERTEX";
    my $aux3 = get_known_entity(8);

    if( scalar @nodes > 0 ) {
      print GRD "1 $::nodes 0 $x $y $z\n";
    }

    $x = get_known_entity(10);
    $y = get_known_entity(20);
    $z = get_known_entity(30);

    if( scalar @nodes == 0 ) {
      $x0 = $x;
      $y0 = $y;
    }

    $::nodes++;
    push(@nodes,$::nodes);
    #print GRD "1 $::nodes 0 $x $y $z\n";
  }

  if( $x == $x0 and $y == $y0 ) {
    pop(@nodes);
    push(@nodes,$nodes[0]);
    #print STDERR "removing last node and closing line...\n";
  } else {
    print GRD "1 $::nodes 0 $x $y $z\n";
  }

  print GRD "\n";

  write_line(\@nodes);
}

sub parse_POINT {

  my $aux1 = get_known_entity(8);

  my $x = get_known_entity(10);
  my $y = get_known_entity(20);
  my $z = get_known_entity(30);

  $::nodes++;
  print GRD "1 $::nodes 0 $x $y $z\n";
}




