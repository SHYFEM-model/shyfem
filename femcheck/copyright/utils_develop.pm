#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub split_developers
{
  my ($name,$year) = @_;

  my @new = ();

  my $ln = length($name);

  my $lymax = 50 - $ln;

  while( length($year) > $lymax ) {
    my $years = "";
    while( length($years) < $lymax ) {
      $year =~ s/([^,]+),//;
      $years .= "$1,";
    }
    $years =~ s/,$//;
    #push(@new,"$name $new");
    push(@new,"$years  $name");
  }
  push(@new,"$year  $name");

  return @new;
}

sub subst_developers
{
  my $rev = shift;

  return $rev unless @$rev;

  foreach my $line (@$rev) {
    my ($date,$devs,$text) = parse_revision_line($line);
    my $subst = 0;
    my @devs = split(/\&/,$devs);
    foreach my $dev (@devs) {
      if( $::subst_dev_names{$dev} ) {
        my $new = $::subst_dev_names{$dev};
	print "   substituting developer $dev with $new\n";
	$dev = $new;
	$subst++;
      }
    }
    if( $subst ) {
      $devs = join("&",@devs);
      $line = "$::comchar $date\t$devs\t$text";
    }
  }

  return $rev;
}

sub count_developers
{
  my $rev = shift;

  make_dev_names();

  return 0 unless @$rev;

  foreach (@$rev) {
    my ($date,$devs,$text) = parse_revision_line($_);
    #print "$date,$devs,$text\n";
    if( $date ) {
      insert_developers($devs,$date);
    }
  }

  my @keys = keys %::devcount;
  my $n = @keys;

  return $n;
}

sub insert_developers
{
  my ($dev,$date) = @_;

  my $year = $date;
  $year =~ s/^.*\.//;

  my @devs = split(/\&/,$dev);
  foreach my $d (@devs) {
    #print STDERR "$d   $date  $::file\n" unless $::copy;
    unless( $::dev_names{$d} ) {
      print STDERR "*** in file $::file no such developer: $d\n";
    }
    $::devyear{$d} .= "$year,";
    $::devcount{$d}++;
  }
}

sub handle_developers
{
  my ($rev,$copy) = @_;

  make_dev_names();

  return 0 unless @$rev;
  return 0 unless @$copy;

  #--------------------------------------------
  # parse developer list
  #--------------------------------------------

  my @keys = keys %::devcount;
  my $n = @keys;
  #print "  $::file: $n\n";
  foreach my $key (@keys) {
    my $years = $::devyear{$key};
    my ($newyears,$nyears,$first) = compact_years($years);
    $::devyear{$key} = $newyears;
    $::devnyear{$key} = $nyears;
    $::devfirst{$key} = $first;
    my $name = $::dev_names{$key};
    unless( $name ) {
      print STDERR "*** no such developer in file $::file: $key\n";
      $name = "unknown";
    }
    $::devname{$key} = $name;
  }

  #--------------------------------------------
  # sort copyright lines and split long lines
  #--------------------------------------------

  my @newcopy = ();
  foreach my $key (sort { 
			   $::devfirst{$a} <=> $::devfirst{$b} 
  			or $::devnyear{$b} <=> $::devnyear{$a} 
  			or $::devname{$b} cmp $::devname{$a} 
			} keys %::devcount ) {
    my $name = $::devname{$key};
    next if $name eq "unknown";
    my $year = $::devyear{$key};
    #my @lines = split_developers("$name ($key)",$year);
    my @lines = split_developers("$name",$year);
    push(@newcopy,@lines);
  }

  #--------------------------------------------
  # write complete copyright lines
  #--------------------------------------------

  foreach (@newcopy) { 
    $_ = $::comchar . "    ". "Copyright (C) $_";
  }

  #--------------------------------------------
  # end of routine
  #--------------------------------------------

  return \@newcopy;
}

sub compact_years
{
  my $years = shift;

  $years =~ s/,$//;			#pop trailing ,
  #print "$years\n";
  my @years = split(",",$years);

  my $last = 0;
  my @new = ();
  foreach (@years) {
    push(@new,$_) if $_ != $last;
    $last = $_;
  }
  #print join(",",@new),"\n";
  my $nyears = @new;
  my $first = $new[0];

  my $line = "";
  my $start = shift(@new);
  my $end = $start;
  foreach (@new) {
    die "internal error... \n" if( $_ == $end );
    if( $_ - $end  == 1 ) {
      $end = $_;
    } else {
      if( $start == $end ) {
        $line .= "$start,";
      } else {
        $line .= "$start-$end,";
      }
      $start = $_;
      $end = $_;
    }
  }
  if( $start == $end ) {
    $line .= "$start,";
  } else {
    $line .= "$start-$end,";
  }

  $line =~ s/,$//;			#pop trailing ,
  #print "$line\n";
  return ($line,$nyears,$first);
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub make_dev_names {

  return if $::dev_initialized;
  $::dev_initialized = 1;

  %::dev_names = (
         'ggu' => 'Georg Umgiesser'
        ,'aac' => 'Andrea Cucco'
        ,'aar' => 'Aaron Roland'
        ,'ccf' => 'Christian Ferrarin'
        ,'cpb' => 'unknown'
        ,'dbf' => 'Debora Bellafiore'
        ,'dmc' => 'Donata Melaku Canu'
        ,'erp' => 'Erik Pascolo'
        ,'fdp' => 'Francesca De Pascalis'
        ,'gir' => 'Ginevra Rosati'
        ,'clc' => 'Celia Laurent'
        ,'cla' => 'Carl Amos'
        ,'isa' => 'Isabella Scroccaro'
        ,'ivn' => 'Ivan Federico'
        ,'laa' => 'Leslie Aveytua'
        ,'lcz' => 'Lucia Zampato'
        ,'mbj' => 'Marco Bajo'
        ,'mcg' => 'Michol Ghezzo'
        ,'pzy' => 'Petras Zemlys'
        ,'unm' => 'Urs Neumeier'
        ,'wmk' => 'William McKiver'
        ,'riz' => 'Rasa Idzelyte'
        ,'git' => 'Git Versioning System'
    );

  %::subst_dev_names = (
         'georg' => 'ggu'
        ,'dmk' => 'dmc'
        ,'cl' => 'clc'
        ,'gr' => 'gir'
    );
}

#--------------------------------------------------------------
1;
#--------------------------------------------------------------

