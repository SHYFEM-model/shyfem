#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# grd utility routines -> old routines , do not use

##############################################################

sub get_nodes {

  $_ = shift;
  s/^\s*//;
  return split;
}

sub put_nodes {

  return join(" ",@_);
}

##############################################################

sub writegrd {

  print "\n";

  foreach $line (@comments) {
    print "$line\n";
  }
  print "\n";
    
  for($i=1;$i<=$nodes;$i++) {
    print "1 $nn[$i] $nt[$i] $nx[$i] $ny[$i] $nh[$i]\n";
  } 
  print "\n";

  for($i=1;$i<=$elems;$i++) {
    print "2 $en[$i] $et[$i] $eg[$i] $ee[$i] $eh[$i]\n";
  } 
  print "\n";

  for($i=1;$i<=$lines;$i++) {
    print "3 $ln[$i] $lt[$i] $lg[$i]\n";
    &write_nodes($le[$i]);
  } 
  print "\n";
}

sub write_nodes {

  my $list = shift;
  my $j = 0;

  my @f = &get_nodes($list);

  foreach (@f) {
    $j++;
    print " $_";
    print "\n" if $j%10 == 0;
  }
  print "\n";

}

###########################################################

sub readgrd {

  my $file = $_[0];

  open(FILE,"$file");

  while( $_ = &nextitem ) {

	@f = split;

	$item = $f[0];

	if( $item == 0 ) {
		&insert_comment;
	} elsif( $item == 1 ) {
		&insert_node;
	} elsif( $item == 2 ) {
		&insert_elem;
	} elsif( $item == 3 ) {
		&insert_line;
	} else {
		die "Unknown item: $_\n";
	}
  }

  close(FILE);

  %nn = &make_hash( @nn );
  %en = &make_hash( @en );
  %ln = &make_hash( @ln );

  &check_hash(\@nn,\%nn);
  &check_hash(\@en,\%en);
  &check_hash(\@ln,\%ln);
}

###################################

sub nextitem {

  my $newline;

  $line = &getline;

  while( $newline = &getline ) {
    #print "$newline\n";
    if( $newline =~ /^\s*$\n/ ) {	#empty line
	last;
    } elsif( $newline =~ /^\s+/ ) {		#conti line
	$line .= " $newline";
    } else {
	&ungetline($newline);		#save for next call
	last;
    }
  }

  $line =~ s/\n/ /g;
  return $line;
}

###################################

sub getline {

  if( @oldlines ) {
    return shift(@oldlines);
  } else {
    return <FILE>;
  }
}

sub ungetline {

  push(@oldlines,$_[0]);
}

###################################

sub insert_comment {
  push(@comments,$_);
}

sub insert_node {
  $nodes++;
  $nn[$nodes] = $f[1];
  $nt[$nodes] = $f[2];
  $nx[$nodes] = $f[3];
  $ny[$nodes] = $f[4];
  $nh[$nodes] = $f[5];
}

sub insert_elem {
  $elems++;
  $en[$elems] = $f[1];
  $et[$elems] = $f[2];
  $eg[$elems] = $f[3];
  my $aux = "";
  my $n = $eg[$elems];
  my $j = 3;
  for($i=0;$i<$n;$i++) {
    $j++;
    $aux .= " $f[$j]";
  }
  $ee[$elems] = $aux;
  $j++;
  $eh[$elems] = $f[$j];
}

sub insert_line {
  $lines++;
  $ln[$lines] = $f[1];
  $lt[$lines] = $f[2];
  $lg[$lines] = $f[3];
  my $aux = "";
  my $n = $lg[$lines];
  my $j = 3;
  for($i=0;$i<$n;$i++) {
    $j++;
    $aux .= " $f[$j]";
  }
  $le[$lines] = $aux;
  $j++;
  $lh[$lines] = $f[$j];
}

###################################

sub make_hash {

  my $n = @_;
  $n--;
  my %hash = {};
  my $i;

  for($i=1;$i<=$n;$i++) {
    my $num = $_[$i];
    $hash{$num} = $i;
  }

  return %hash;
}

sub check_hash {

  my $ra = shift;
  my $rh = shift;

  my $n = @$ra;
  $n--;
  my $i;
  my $num;

  for($i=1;$i<=$n;$i++) {
    $num = $$ra[$i];
    if( $$rh{$num} != $i ) {
      die "Inconsistent hash (1): $i $num $$rh{$num}\n";
    }
  }

  my @hk = keys(%$rh);

  foreach $num (@hk) {
    my $ind = $$rh{$num};
    if( $$ra[$ind] != $num ) {
      die "Inconsistent hash (2): $ind $num $$ra{$ind}\n";
    }
  }
}

###################################
1;
###################################

