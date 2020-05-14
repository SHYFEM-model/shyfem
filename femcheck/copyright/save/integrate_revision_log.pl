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
$::obsolete = 0;		#check for obsolete date

$::extract = 0 unless $::extract;

my $in_revision = 0;
my $has_revision_log = 0;

$::names = ();
make_names();

my $file = $ARGV[0];
my $revlog = ReadRevisionLog("revlog.tmp");
my $revlog_git = ReadRevisionLog("revlog.git.tmp");

if( not $revlog and not $revlog_git ) {
  die "no revision log found at all in $file\n";
} elsif( not $revlog_git ) {
  CopyFile();
} elsif( not $revlog ) {
  InsertRevisionLog($revlog_git);
} else {
  my $revlog_new = CombineRevisionLog($revlog,$revlog_git);
  SubstituteRevisionLog($revlog_new);
}

#---------------------------------------------------------

sub ReadRevisionLog
{
  my $file = shift;

  my @items = ();

  if( not open(REV,"<$file") {
    #print STDERR "No revision log in file $file\n";
    return;
  }

  while(<REV>) {
    chomp;
    next if /^\s*$/;		#skip empty line
    next if /^[cC!*]\s*$/;	#skip just comment line
    s/\^[cC!*]\s*//;		#strip comment in front
    my %item = ();
    if( /^(\S+)\s+(\S+)\s+(.*)/ ) {
      $item{date} = $1;
      $item{name} = $2;
      $item{text} = $3;
      $item{idate} = make_idate($1);
      push(@items,\%item);
    } else {
      print STDERR "cannot parse: $_\n";
    }
  }

  close(REV);

  return \@items;
}

#---------------------------------------------------------

sub CopyFile
{
  while(<>) {
    print;
  }
}

#---------------------------------------------------------

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

sub make_idate {

  my $date = shift;

  if( $date =~ /^(\d\d)\.(\d\d)\.(\d\d\d\d)$/ ) {
    my $idate = $1 + 100*$2 + 10000*$3;
    return $idate
  } else {
    die "cannot parse date: $date\n";
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

#---------------------------------------------------------

