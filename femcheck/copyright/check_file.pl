#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femcheck/copyright"
	,"$ENV{HOME}/shyfem/femcheck/copyright");

use strict;
use revision_log;

my $file = $ARGV[0];

$::warn = 0 unless $::warn;		#warn for revision log out of section
$::obsolete = 0 unless $::obsolete;	#check for obsolete date
$::extract = 0 unless $::extract;	#extract revision log
$::integrate = "" unless $::integrate;	#integrate revision log

#------------------------------------------------------------------------

init_revision();

#print STDERR "integrate: $::integrate\n";

if( $::integrate ) {
  integrate_revision_log($file,$::integrate);
  exit 1;
}

my $revlog = get_revision_log();

check_dates_of_revision_log($revlog,$file);

if( $::extract and $::has_revision_log ) {
  write_revision_log("revlog.tmp",$revlog);
}

check_file();

my $return = $::has_revision_log;
$return = $::is_manual if $::is_manual;

exit $return;

#------------------------------------------------------------------------

sub check_file
{
  if( $::has_copyright == 0 ) {
    print STDERR "file has no copyright: $file\n";
  }

  if( $::has_shyfem == 0 ) {
    print STDERR "file has no shyfem line: $file\n";
  }

  if( $::is_manual == 1 ) {
    print STDERR "file has manual copyright: $file\n";
  }
}

#------------------------------------------------------------------------

