#!/usr/bin/perl -s -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# parses STR file and writes info to stdout or grd file
#
# possible command line options: see subroutine FullUsage
#
#--------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");

use str;
use grd;
use strict;

#-------------------------------------------------------------
# command line options
#-------------------------------------------------------------
$::h = 0 unless $::h;
$::help = 0 unless $::help;
$::quiet = 0 unless $::quiet;
$::bnd = 0 unless $::bnd;
$::files = 0 unless $::files;
$::zip = 0 unless $::zip;
$::rewrite = 0 unless $::rewrite;
$::sect = "" unless $::sect;
$::txt = 0 unless $::txt;
$::debug = 0 unless $::debug;
$::value = "" unless $::value;
$::replace = "" unless $::replace;
#-------------------------------------------------------------

#-------------------------------------------------------------
# info on names in sections
#-------------------------------------------------------------
@::bound_names = qw/ boundn conzn tempn saltn
		bio2dn sed2dn tox3dn
		bfm1bc bfm2bc bfm3bc /;
@::name_names = qw/ bound wind rain qflux ice restrt gotmpa
		bio bios toxi
		conzin saltin tempin /;
@::aquabc_names = qw/ biocon bioscon biolight bioaow 
		bioaos bioph biotemp bioload /;
@::lagrg_names = qw/ lagra lagrf lagrt /;
@::sedtr_names = qw/ sedp sedt sedcon /;
#-------------------------------------------------------------

my $file = $ARGV[0];

my $str = new str;
$str->{quiet} = $::quiet;
$str->read_str($file) if $file;

if( $::h or $::help ) {
  FullUsage();
} elsif( not $file ) {
  Usage();
} elsif( $::bnd ) {
  show_bnd_nodes($str);
} elsif( $::files ) {
  show_files($str);
} elsif( $::zip ) {
  my $files = show_files($str);
  push(@$files,$file);		#add str-file to archive
  zip_files($files);
} elsif( $::rewrite ) {
  #$str->print_sections();;
  $str->write_str("new.str");;
  print STDERR "str-file written to new.str\n";
} elsif( $::sect ) {
  show_sect($str,$::sect);
} elsif( $::value ) {
  if( $::replace ) {
    replace_value($str,$::value,$::replace);
    $str->write_str("replace.str");;
    print STDERR "str-file written to replace.str\n";
  } else {
    my $val = show_value($str,$::value);
    $val = "(unknown)" unless $val;
    #print "$::value = $val\n";
    print "$val\n";
  }
} else {
  Usage();
}


#------------------------------------------------------------

sub Usage {

  print STDERR "Usage: strparse.pl [-h|-help] [-options] str-file\n";
  exit 0;
}

sub FullUsage {

  print STDERR "Usage: strparse.pl [-h|-help] [-options] str-file\n";
  print STDERR "  options:\n";
  print STDERR "    -h!-help      this help screen\n";
  print STDERR "    -quiet        be as quiet as possible\n";
  print STDERR "    -bnd          extract boundary nodes\n";
  print STDERR "    -files        extract names of forcing files\n";
  print STDERR "    -zip          zips forcing files, grid, str in one file\n";
  print STDERR "    -rewrite      rewrite the str file\n";
  print STDERR "    -value=var    show value of var ([sect:]var)\n";
  print STDERR "    -replace=val  replace value of var with val and rewrite\n";
  #print STDERR "    -sect=sect    writes contents of section\n";
  print STDERR "    -txt          write nodes as text and not in grd format\n";
  print STDERR "  if -replace is given also -value must be specified\n";
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

  my $outfile;
  if( $::txt ) {
    $outfile = "bnd_str.txt";
  } else {
    $outfile = "bnd_str.grd";
  }
  open(OUT,">$outfile");

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "bound" ) {
      show_nodes($str,$grid,$sect);
    }
  }

  close(OUT);
  print STDERR "nodes written to file $outfile\n";
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
    if( ref($value) eq "ARRAY" ) {
      write_array_new($value,"$sect_id :  kbound = \n");
      @list = @$value;
    } else {
      print "$sect_id :    kbound = $value\n";
      push(@list,$value);
    }
  }
  
  if( $::txt ) {
    write_nodes_to_txt($sect_number,\@list);
  } else {
    #write_nodes_to_grd($ibtyp,$grid,\@list);
    write_nodes_to_grd($sect_number,$grid,\@list);
  }
}

sub write_nodes_to_txt {

  my ($id,$list) = @_;

  my $n = @$list;
  print OUT "$id  $n\n";
  while( $n-- ) {
    my $val = shift(@$list);
    print OUT "$val\n";
  }
}

sub write_nodes_to_grd {

  my ($ibtyp,$grid,$list) = @_;

  return if( $ibtyp <= 0 );
  my $type = $ibtyp;

  foreach my $node (@$list) {
    my $item = $grid->get_node($node);
    my $x = $item->{x};
    my $y = $item->{y};
    print OUT "1 $node $type $x $y\n";
  }
}

sub show_files {

  my $str = shift;

  my @files = ();
  my $new = "";

  my $basin = $str->get_basin();
  $basin .= ".grd";
  $basin =~ s/^\s+//;
  print "basin :    grid = '$basin'\n";
  push(@files,$basin);

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};
    my $section_id = $sect->{id};

    $new = "";
    if( $sect->{name} eq "bound" ) {
      $new = show_name($str,$sect,\@::bound_names);
    } elsif( $sect->{name} eq "name" ) {
      $new = show_name($str,$sect,\@::name_names);
    } elsif( $sect->{name} eq "aquabc" ) {
      $new = show_name($str,$sect,\@::aquabc_names);
    } elsif( $sect->{name} eq "lagrg" ) {
      $new = show_name($str,$sect,\@::lagrg_names);
    } elsif( $sect->{name} eq "sedtr" ) {
      $new = show_name($str,$sect,\@::sedtr_names);
    }

    push(@files,@$new) if $new;
  }

  return \@files;
}

sub replace_value {

  my ($str,$name,$replace,$sect_name,$sect_number) = @_;

  #$sect_name = "para" unless $sect_name;
  #$sect_number = "" unless $sect_number;

  $str->set_value($name,$replace,$sect_name,$sect_number);
}

sub show_value {

  my ($str,$name,$sect_name,$sect_number) = @_;

  #$sect_name = "para" unless $sect_name;
  #$sect_number = "" unless $sect_number;

  return $str->get_value($name,$sect_name,$sect_number);
}

sub show_name {

  my ($str,$sect,$var_names) = @_;

  my @files = ();

  my $sect_name = $sect->{name}; 
  my $sect_number = $sect->{number}; 
  my $sect_id = $sect->{id}; 

  if( $::debug ) {
    print STDERR "searching files: $sect_name ($sect_number) $sect_id\n";
  }

  foreach my $name (@$var_names) {
    my $value = $str->get_value($name,$sect_name,$sect_number);
    print STDERR "  $name  (no value)\n" if not defined $value and $::debug;
    if( defined $value ) {
      print STDERR "  $name  $value\n" if $::debug;
      print "$sect_id :    $name = $value\n";
      $value =~ s/\'//g;
      push(@files,$value);
    }
  }

  return \@files;
}

sub write_array_new {

  my ($array,$extra) = @_;

  my $nval = 10;		# how many values on one line

  print "$extra" if $extra;

  my $i = 0;
  foreach my $item (@$array) {
    $i++;
    print " $item";
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

sub show_sect {

  my ($str,$sname) = @_;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    my $nnn = $sect->{name}; print "section: $nnn\n";

    if( $sect->{name} eq "$sname" ) {
        #show_name($str,$sect,\@::bound_names);
	my $array = $sect->{array};
	foreach my $numb (@$array) {
	  print "$numb\n";
	}
    }
  }
}

sub zip_files {

  my ($files) = @_;

  my $zipfile = "new_str_zip.zip";

  unlink "$zipfile" if ( -f $zipfile );

  foreach my $file (@$files) {
    print STDERR "zip file: $file\n";
    system("zip $zipfile $file");
  }

  print STDERR "files have been zipped into file $zipfile\n";
}

