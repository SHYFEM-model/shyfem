#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

%::names = ();


make_names();

while(<>) {

  chomp;
  check_old_revision();
  print "$_\n";
}

#--------------------------------------------------------------

sub check_old_revision {

  if( /^[cC!].*(\d{2}\.\d{2}\.\d{2})\s+/ ) {
    my $orig = $_;
    if( /^([cC!])\s+(.+)\s+(\d{2}\.\d{2}\.\d{2})\s+(.+)/ ) {
      my $comment = $1;
      my $verb = $2;
      my $date = $3;
      my $descrp = $4;
      $verb =~ s/\s+on$//;
      $date = adjust_year($date);
     
      if( $verb eq "written" or $verb eq "revised" or $verb eq "changed" ) {
        ;
      } else {
	print STDERR "*** cannot parse (1): $_\n";
      }
      $descrp =~ s/\s*by ggu\s*//;
      print STDERR "$orig\n";
      $_ = "$comment $date ggu\t$descrp";
      print STDERR "$_\n";
    } else {
      print STDERR "*** cannot parse (2): $_\n";
    }
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

