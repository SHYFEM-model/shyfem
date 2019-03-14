#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use strict;

$::count = 0 unless $::count;
$::copy = 0 unless $::copy;
$::subst = 0 unless $::subst;
$::old = 0 unless $::old;
$::new = 0 unless $::new;
$::revision = 0 unless $::revision;	#exit if no revision log
$::file = $ARGV[0];
%::devyear = ();
%::devcount = ();
%::names = ();
$::file = $ARGV[0];

$::write = 1;
$::write = 0 if $::copy;
$::write = 0 if $::subst;

$::copystart = '!    Copyright (C)';
$::copystd   = '!    Copyright (C) 1985-2019  Georg Umgiesser';
$::copyold   = '';
$::copynew   = '';
$::copydone  = 0;

make_names();

count() if $::count;	# this just counts developers and writes it to terminal

my $newcopy = read_new_copyright() if $::subst;

while(<>) {

  chomp;
  check_copyright() if $::copy;
  subst_copyright($newcopy) if $::subst;

  if( $::old ) {
    check_old_revision();
  } else {
    check_new_revision();
  }
}

if( $::copy ) {
  make_copyright();
  handle_copyright();
}

#--------------------------------------------------------------

sub handle_copyright {

  if( $::copyold ne $::copynew ) {
    print "----------- old copyright -------------- $::file \n";
    print "$::copyold";
    print "----------- new copyright --------------\n";
    print "$::copynew";
    print "----------- end copyright --------------\n";
  }
}

sub subst_copyright
{

  my $newcopy = shift;

  if( /Copyright \(C\)/ ) {
    return if $::copydone;
    $::copydone = 1;
    print $newcopy ;
  } else {
    print "$_\n";
  }
}

sub check_copyright
{
  if( /Copyright \(C\)/ ) {
    $::copyold .= "$_\n";
  }
}

sub make_copyright
{
  my %hash = %::devcount;
  my @keys = sort_copyright(\%hash);
  my $n = 0;

  foreach my $key (@keys) {
     my $count = $::devcount{$key};
     my $years = consolidate_years($::devyear{$key});
     my $name = $::names{$key};
     unless( $name ) {
	print STDERR "*** no such name: $key\n";
	$name = "unknown";
     }
     print_copyright($years,$name);
     #print "$key  $count   $years\n"; 
     $n++;
  }

  if( $n == 0 ) {
    print STDERR "*** no revision log in $::file\n";
    $::copynew = "$::copystd\n";
  }

  exit(0) if( $::revision );
}

sub sort_copyright {

  my $rhash = shift;

  my @keys = sort { $rhash->{$b} <=> $rhash->{$a} } keys %$rhash;

  my $n = @keys;
  return () unless $n;

  my $changed = 1;
  while( $changed ) {
    my @old = @keys;
    @keys = ();
    $changed = 0;
    my $a = $old[0];
    for( my $i=1; $i<$n; $i++ ) {
      my $b = $old[$i];
      my $ayear = $::devyear{$a};
      $ayear =~ s/^(\d+).*/$1/;
      my $byear = $::devyear{$b};
      $byear =~ s/^(\d+).*/$1/;
      if( $ayear > $byear ) {
        push(@keys,$b);
      } elsif( $ayear == $byear ) {
	if( $::names{$a} gt $::names{$b} ) {
          push(@keys,$b);
	}
      } else {
        push(@keys,$a);
        $a = $b;
      }
    }
    push(@keys,$a);
  }
  return @keys;
}

sub print_copyright {

  my ($years,$name) = @_;

  my $max = 75;
  my $ly = length($years);
  my $ln = length($name);
  my $lc = length($::copystart);

  if( $lc+$ly+$ln > $max ) {
    print STDERR "*** must break down...: $years\n";
    $max = $max - $lc - $ln;
    my @y = split(/,/,$years);
    my $line = "";
    foreach my $y (@y) {
      my $n = length($y);
      if( length($line) + $n + 1 > $max ) {
        $::copynew .= "$::copystart $line  $name\n"; 
        $line = "";
      }
      $line .= "," if $line;
      $line .= $y;
    }
    $::copynew .= "$::copystart $line  $name\n"; 
  } else {
    $::copynew .= "$::copystart $years  $name\n"; 
  }
}

sub consolidate_years {

  my $years = shift;

  $years =~ s/,$//;
  my @years = split(/,/,$years);
  my @new = ();

  my $oldyear = shift(@years);
  push(@new,$oldyear);
  foreach my $year (@years) {
    next if $year == $oldyear;
    push(@new,$year);
    $oldyear = $year;
  }

  my $line = join(",",@new);
  #print "$line\n";

  $line = "";
  my $start = shift(@new);
  my $end = $start;
  foreach my $year (@years) {
    if( $year - $end > 1 ) {
      if( $start == $end ) {
        $line .= "$end,";
      } else {
        $line .= "$start-$end,";
      }
      $start = $year;
      $end = $year;
    } else {
      $end = $year;
    }
  }
  if( $start == $end ) {
    $line .= "$end,";
  } else {
    $line .= "$start-$end,";
  }
  $line =~ s/,$//;

  return $line;
}

#--------------------------------------------------------------

sub check_old_revision {

  if( /^[cC!].*(\d{2}\.\d{2}\.\d{2})\s+/ ) {
    print "$_\n";
  }
}

sub check_new_revision {

  if( /^[cC!]\s+(\d{2}\.\d{2}\.\d{4})\s+(\S+)\s+/ ) {
    my $dev = $2;
    my $date = $1;
    my $year = $date;
    $year =~ s/^.*\.//;
    my @devs = split(/\&/,$dev);
    foreach my $d (@devs) {
      print "$d   $date  $::file\n" if $::write;
      $::devyear{$d} .= "$year,";
      $::devcount{$d}++;
    }
  }
}

#--------------------------------------------------------------

sub read_new_copyright {

  my $file = "copy.tmp";
  my $in_new = 0;
  my $line = "";

  open(FILE,"<$file") || die "Cannot open file: $file\n";
  while(<FILE>) {
    $in_new = 0 if( /-- end copyright --/ );
    $in_new = 1 if( /-- new copyright --/ );
    next if( /^---/ );
    $line .= $_ if( $in_new );
  }
  close(FILE);

  #print STDERR "copyright read:\n$line";

  return $line;
}

#--------------------------------------------------------------

sub count
{
  my %dev = ();

  while(<>) {
    my @f = split;
    my $d = $f[0];
    $dev{$d}++;
  }

  foreach my $d (sort keys %dev) {
    my $c = $dev{$d};
    my $name = $::names{$d};
    print "$d  $c  $name\n";
  }

  exit 0;
}

#--------------------------------------------------------------

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

