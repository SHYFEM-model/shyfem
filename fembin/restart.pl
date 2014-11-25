#!/usr/bin/perl -s -w
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

my $file = $ARGV[0];

my $str = new str;
$str->read_str($file) if $file;

if( $::h or $::help ) {
  FullUsage();
} elsif( not $file ) {
  Usage();
} else {
  my $new_str = prepare_restart($str,$file);
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
  #print STDERR "    -sect=sect    writes contents of section\n";
  exit 0;
}

#------------------------------------------------------------

sub prepare_restart {

  my ($str,$file) = @_;

  my $simul = $str->get_simul();
  $file =~ s/\.str$//;

  my ($new_simul,$new_str,$post) = make_new_filenames($simul,$file);

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
  print STDERR "please restart simulation as \"ht < $new_str\"\n";

  return $new_str;
}

sub make_new_filenames {

  my ($simul,$str) = @_;

  my $new_simul;
  my $new_str;
  my $post;
  my $ok = 0;

  foreach my $i (1..99) {
    $post = "_rst_$i";
    $new_simul = $simul . $post;
    $new_str = $str . $post . ".str";
    next if ( -f $new_str );
    next if ( -f "$simul.inf" );
    next if ( -f "$simul.ous" );
    next if ( -f "$simul.nos" );
    $ok = 1;
    last;
  }

  die "Cannot find names to use for restart... sorry.\n" unless $ok;

  return ($new_simul,$new_str,$post);
}

