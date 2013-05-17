#!/usr/bin/perl -s -w
#
# parses STR file and writes info to stdout or grd file
#
# possible command line options:
#
#	-bnd	open boundary nodes (writes bnd_str.grd)
#	-files	files with open boundary conditions
#	-txt	writes boundary nodes in txt format
#
#--------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use str;
use grd;
use strict;

#-------------------------------------------------------------
# command line options
#-------------------------------------------------------------
$::files = 1 if $::files;
$::bnd = 1 if $::bnd;
$::txt = 1 if $::txt;
#-------------------------------------------------------------

#-------------------------------------------------------------
@::bound_names = qw/ boundn conzn tempn saltn
		bio2dn sed2dn tox3dn
		bfm1bc bfm2bc bfm3bc /;
@::name_names = qw/ bound wind rain qflux restrt gotmpa
		bio bios toxi
		conzin saltin tempin /;
@::aquabc_names = qw/ biocon bioscon biolight bioaow 
		bioaos bioph biotemp bioload /;
@::lagrg_names = qw/ lagra lagrf lagrt /;
@::sedtr_names = qw/ sedp sedt sedcon /;
#-------------------------------------------------------------

my $file = $ARGV[0];
Usage() unless $file;

my $str = new str;
$str->read_str($file);

if( $::bnd ) {
  show_bnd_nodes($str);
} elsif( $::files ) {
  show_files($str);
} else {
  Usage();
}

#------------------------------------------------------------

sub Usage {

  print STDERR "Usage: strparse.pl {-bnd|-files} [-txt] str-file\n";
  exit 0;
}

#------------------------------------------------------------

sub show_bnd_nodes {

  my $str = shift;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  my $basin = $str->get_basin();
  my $grid = new grd;
  $grid->readgrd("$basin.grd");

  open(GRD,">bnd_str.grd");

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "bound" ) {
      show_nodes($str,$grid,$sect);
    }
  }

  close(GRD);
  print STDERR "nodes written to file bnd_str.grd\n";
}

sub show_nodes {

  my ($str,$grid,$sect) = @_;

  my $sect_name = $sect->{name}; 
  my $sect_number = $sect->{number}; 
  my $sect_id = $sect->{id}; 

  my @list = ();

  my $ibtyp = $str->get_value('ibtyp',$sect_name,$sect_number);
  $ibtyp = 1 if not defined $ibtyp;
  my $value = $str->get_value('kbound',$sect_name,$sect_number);
  if( defined $value ) {
    if( $::txt ) {
      write_txt($sect_number,$value);
    } elsif( ref($value) eq "ARRAY" ) {
      write_array_new($value,"$sect_id :  kbound = \n");
      @list = @$value;
    } else {
      print "$sect_id :    kbound = $value\n";
      push(@list,$value);
    }
  }
  
  write_nodes_to_grd($ibtyp,$grid,\@list);
}

sub write_nodes_to_grd {

  my ($ibtyp,$grid,$list) = @_;

  return if( $ibtyp <= 0 );
  my $type = $ibtyp;

  foreach my $node (@$list) {
    my $item = $grid->get_node($node);
    my $x = $item->{x};
    my $y = $item->{y};
    print GRD "1 $node $type $x $y\n";
  }
}

sub show_files {

  my $str = shift;

  my $basin = $str->get_basin();
  print "basin : $basin\n";

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "bound" ) {
      show_name($str,$sect,\@::bound_names);
    } elsif( $sect->{name} eq "name" ) {
      show_name($str,$sect,\@::name_names);
    } elsif( $sect->{name} eq "aquabc" ) {
      show_name($str,$sect,\@::aquabc_names);
    } elsif( $sect->{name} eq "lagrg" ) {
      show_name($str,$sect,\@::lagrg_names);
    } elsif( $sect->{name} eq "sedtr" ) {
      show_name($str,$sect,\@::sedtr_names);
    }

  }

}

sub show_name {

  my ($str,$sect,$var_names) = @_;

  my $sect_name = $sect->{name}; 
  my $sect_number = $sect->{number}; 
  my $sect_id = $sect->{id}; 

  foreach my $name (@$var_names) {
    my $value = $str->get_value($name,$sect_name,$sect_number);
    if( defined $value ) {
      print "$sect_id :    $name = $value\n";
    }
  }
}

sub write_array_new {

  my ($array,$extra) = @_;

  my $nval = 10;		# how many values on one line

  print "$extra" if $extra;

  my $i = 0;
  foreach my $item (@$array) {
    $i++;
    print "   $item";
    print "\n" if $i%$nval == 0;
  }
  print "\n" unless $i%$nval == 0;
}

sub write_txt {

  my ($id,$value) = @_;

  if( ref($value) ne "ARRAY" ) {
    my @value = ();
    $value[0] = $value;
    $value = \@value;
  }

  my $n = @$value;
  print "$id  $n\n";
  while( $n-- ) {
    my $val = shift(@$value);
    print "$val\n";
  }
}

