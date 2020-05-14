#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

$::warn = 0;			#warn for revision out of revision log section
$::obsolete = 0;		#check for obsolete date in text

$::extract = 0 unless $::extract;	#extract revlog to revlog.tmp

my $in_revision = 0;
my $has_revision_log = 0;

$::names = ();
make_names();

open(REV,">revlog.tmp") if $::extract;

while(<>) {

  chomp;

  if( $in_revision == 0 ) {
    if( /revision log :/) {		#start of revision log
      $in_revision = 1;
      $has_revision_log++;
    } else {
      check_revision(0);
    }
  } else {
    if( /^\s*$/ ) {			#empty line - end of revision log
      $in_revision = 0;
    } elsif( /^[!cC]\*\*\*/ ) {		#c*** line
      $in_revision = 0;
    } elsif( /^[!cC]\=\=\=/ ) {		#c=== line
      $in_revision = 0;
    } elsif( /^[!cC].*\s:\s*$/ ) {	#other comment line
      $in_revision = 0;
    } else {
      check_revision(1);
    }
    if( $in_revision ) {		#write rev log to extra file
      print REV "$_\n" if $::extract;
    }
  }

  print "$_\n";
}

close(REV) if $::extract;

exit $has_revision_log;

#--------------------------------------------------------------

sub check_revision {

  my $irv = shift;	#0 if not in rev section, 1 if inside

  return if /^[!cC]\s*$/;
  return if /^[!cC]\*\*\*/;

  if( check_new_revision() ) {
    if( $irv == 0 and $::warn ) {
      print STDERR "new revision out of revision log: $_\n";
    }
  } elsif( check_old_revision() ) {
    if( $irv == 0 and $::warn ) {
      print STDERR "old revision out of revision log: $_\n";
    }
  } elsif( $::obsolete and check_obsolete_revision() ) {
    ;
  } else {
    if( $irv ) {
      print STDERR "*** Cannot parse revision log: $_\n";
    }
  }
}

sub check_new_revision {

  if( /^[cC!]\s+(\d{2}\.\d{2}\.\d{4})\s+(\S+)\s+/ ) {
    my $date = $1;
    my $dev = $2;
    handle_developers($dev,$date);
    return 1;
  } elsif( /^[cC!]\s+\.\.\.\s+/ ) {	#continuation line
    return 1;
  } else {
    return 0;
  }
}

sub check_old_revision {

  if( /^[cC!].*(\d{2}\.\d{2}\.\d{2})\s+/ ) {
    my $orig = $_;
    if( /^([cC!])\s+(.+)\s+(\d{2}\.\d{2}\.\d{2})\s+(.+)/ or
	/^([cC!])(\s+)(\d{2}\.\d{2}\.\d{2})\s+(.+)/ ) {
      my $comment = $1;
      my $verb = $2;
      my $date = $3;
      my $descrp = $4;
      $verb =~ s/\s+on\s*$//;
      $date = adjust_year($date);
     
      if( $verb eq "written" or $verb eq "revised" or $verb eq "changed" ) {
        ;
      } elsif( $verb =~ /\s*/ ) {	#simply old date
        ;
      } else {
	print STDERR "*** cannot parse (1): $_\n";
	return 0;
      }
      $descrp =~ s/\s*by ggu\s*//;
      $descrp =~ s/\s*ggu\s*//;
      $_ = "$comment $date\tggu\t$descrp";
      print STDERR "orig: $orig\n";
      print STDERR "new:  $_\n";
      return 1;
    } else {
      print STDERR "*** cannot parse (2): $_\n";
      return 0;
    }
  } else {
    return 0;
  }
}

sub check_obsolete_revision {

  if( /^[cC!].*(\d{1,2}\s+\d{1,2}\s+\d{2,4})/ ) {	# 12 4 2018
    print STDERR "*** obsolete date (1): $_\n";
    return 1;
  } elsif( /^[cC!].*(\d{1,2}\/\d{1,2}\/\d{2,4})/ ) {	# 12/04/2018
    print STDERR "*** obsolete date (2): $_\n";
    return 1;
  } elsif( /^[cC!].*(\d{1,2}-\d{1,2}-\d{2,4})/ ) {	# 12-04-2018
    print STDERR "*** obsolete date (3): $_\n";
    return 1;
  } elsif( /^[cC!].*(\d{1,4}[-\/]\d{1,4}[-\/]\d{1,4})/ ) {# any with - or /
    print STDERR "*** obsolete date (4): $_\n";
    return 1;
  } elsif( /^[cC!].*(\d{1,4}\s+\d{1,4}\s+\d{1,4})/ ) {	# any with spaces
    print STDERR "*** obsolete date (5): $_\n";
    return 1;
  } else {
    return 0;
  }
}

sub adjust_year {

  my $date = shift;

  my @f = split(/\./,$date);
  if( $f[2] > 70 ) {
    $f[2] += 1900;
  } else {
    $f[2] += 2000;
  }

  $date = join(".",@f);
  return $date;
}

sub handle_developers
{
  my ($dev,$date) = @_;

  my $year = $date;
  $year =~ s/^.*\.//;

  my @devs = split(/\&/,$dev);
  foreach my $d (@devs) {
    #print STDERR "$d   $date  $::file\n" unless $::copy;
    #$::devyear{$d} .= "$year,";
    #$::devcount{$d}++;
  }
}

sub make_names {

  %::names = (
         'ggu' => 'Georg Umgiesser'
        ,'aac' => 'Andrea Cucco'
        ,'aar' => 'Aaron Roland'
        ,'ccf' => 'Christian Ferrarin'
        ,'cpb' => ''
        ,'dbf' => 'Debora Bellafiore'
        ,'dmk' => 'Donata Melaku Canu'
        ,'erp' => 'Erik Pascolo'
        ,'fdp' => 'Francesca De Pascalis'
        ,'isa' => 'Isabella Scroccaro'
        ,'ivn' => 'Ivan Federico'
        ,'laa' => 'Leslie Aveytua'
        ,'lcz' => 'Lucia Zampato'
        ,'mbj' => 'Marco Bajo'
        ,'mcg' => 'Michol Ghezzo'
        ,'pzy' => 'Petras Zemlys'
        ,'wmk' => 'William McKiver'
    );
}

