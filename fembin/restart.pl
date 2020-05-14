#!/usr/bin/perl -s -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# parses STR file and re-writes it for restart
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
#-------------------------------------------------------------

my $strfile = $ARGV[0];

$::start = 1 unless $::start;

my $str = new str;
$str->read_str($strfile) if $strfile;

if( $::h or $::help ) {
  FullUsage();
} elsif( not $strfile ) {
  Usage();
} else {
  my $new_str = prepare_restart($str,$strfile);
  $str->write_str($new_str);
}

#------------------------------------------------------------

sub Usage {

  print STDERR "Usage: restart.pl [-h|-help] [-options] str-file\n";
  exit 0;
}

sub FullUsage {

  print STDERR "Usage: restart.pl [-h|-help] [-options] str-file\n";
  print STDERR "  options:\n";
  print STDERR "    -h!-help      this help screen\n";
  print STDERR "    -start=n      start with n as number to append\n";
  #print STDERR "    -sect=sect    writes contents of section\n";
  exit 0;
}

#------------------------------------------------------------

sub prepare_restart {

  my ($str,$strfile) = @_;

  my $simul = $str->get_simul();
  $strfile =~ s/\.str$//;

  my ($new_simul,$new_str,$post) = make_new_filenames($simul,$strfile);

  my $restart_file = $simul . ".rst";

  my $title = $str->get_title();
  $title .= " $post";
  $str->set_title($title);
  $str->set_simul($new_simul);

  $str->set_value("restrt","\'$restart_file\'","name");
  $str->set_value("itrst",-1,"para");

  print STDERR "new str-file name: $new_str\n";
  print STDERR "new simulation name: $new_simul\n";
  print STDERR "restarting from file: $restart_file\n";
  print STDERR "please restart simulation as \"shyfem $new_str\"\n";

  return $new_str;
}

sub make_new_filenames {

  my ($simul,$strfile) = @_;

  my $new_simul;
  my $new_strfile;
  my $post;
  my $ok = 0;
  my $start = $::start;

  if( $simul =~ /^(.+)_rst_(\d+)$/ ) {
    $simul = $1;
    $start = $2 if $2 > $start;
  }

  if( $strfile =~ /^(.+)_rst_(\d+)$/ ) {
    $strfile = $1;
    $start = $2 if $2 > $start;
  }

  foreach my $i ($start..99) {
    #$i = "0" . $i if $i < 10;
    $post = "_rst_$i";
    $new_simul = $simul . $post;
    $new_strfile = $strfile . $post . ".str";
    next if ( -f "$new_strfile" );
    next if ( -f "$new_simul.inf" );
    next if ( -f "$new_simul.ous" );
    next if ( -f "$new_simul.nos" );
    $ok = 1;
    last;
  }

  die "Cannot find names to use for restart... sorry.\n" unless $ok;

  return ($new_simul,$new_strfile,$post);
}

