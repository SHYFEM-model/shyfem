#!/usr/bin/perl -ws
#
# checks include and computes dipendencies for fortran files
#
#------------------------------------------------------

use strict;

$::nomake = 1 if $::nomake;
$::debug = 1 if $::debug;

$::make = 1 unless $::nomake;

my $mkdp = "# DO NOT DELETE THIS LINE -- make depend depends on it.";

my @lines = ();

foreach my $file (@ARGV) {

  print STDERR "$file\n" if $::debug;

  open(FILE,"$file") || die "Cannot open file $file\n";
  my $rlist = [];

  if( $file =~ /\.f$/ ) {
    check_f_file($rlist);
  } elsif( $file =~ /\.c$/ ) {
    check_c_file($rlist);
  } else {
    die "Cannot handle type of file: $file\n";
  }


  close(FILE);

  my $line = write_inc($file,$rlist);
  print "$line\n" if $line;
  push(@lines,$line) if $line;
}

if( not $::make ) {
  print STDERR "Makefile not changed...\n";
  exit 0;
}

print STDERR "changing Makefile...\n";
my $mfile = change_makefile(\@lines);

rename("$mfile","$mfile.bak");
rename("$mfile.new","$mfile");

#----------------------------------------------------------

sub check_f_file {

  my $rlist = shift;

  while( <FILE> ) {
    if( /^\s+include\s*['"]\s*([\w.]+)\s*['"]\s*$/) {
      my $hfile = $1;
      print STDERR "include found: $hfile\n" if $::debug;
      insert_inc($rlist,$hfile);
    }
  }
}

sub check_c_file {

  my $rlist = shift;

  while( <FILE> ) {
    if( /^\s*\#\s*include\s*['"]\s*([\w.]+)\s*['"]\s*$/) {
      my $hfile = $1;
      print STDERR "include found: $hfile\n" if $::debug;
      insert_inc($rlist,$hfile);
    }
  }
}

#----------------------------------------------------------

sub change_makefile {

  my $ra = shift;

  my $mkdp_found = 0;
  my $mfile = "";

  unless( $mfile ) {
    open(MAKE,"<Makefile") && ($mfile = "Makefile");
  }
  unless( $mfile ) {
    open(MAKE,"<makefile") && ($mfile = "makefile");
  }
  unless( $mfile ) {
    die "no Makefile found\n";
  }

  open(NEW,">$mfile.new") || die "cannot open file: $mfile.new\n";

  while(<MAKE>) {
    if( /^$mkdp/ ) {
        print STDERR "make line found...\n";
	$mkdp_found = 1;
        last;
    }
    print NEW;
  }

  print NEW "\n" unless $mkdp_found;		# just for beauty
  print NEW "$mkdp\n";
  print NEW "\n";

  $ra = fit_line($ra);

  foreach my $line (@$ra) {
    print NEW "$line\n";
  }
  print NEW "\n";
  
  close(MAKE);
  close(NEW);

  return $mfile;
}

sub write_inc {

  my ($file,$rlist) = @_;
  
  return unless scalar @$rlist;

  $file =~ s/\.[fc]$/.o:/;
  foreach my $f (@$rlist) {
    $file .= " $f";
  }

  return $file;
}

sub insert_inc {

  my ($rlist,$hfile) = @_;

  foreach my $f (@$rlist) {
    return if $f eq $hfile;		# already inserted
  }

  push(@$rlist,$hfile);
}

sub strip {

  my $s = shift;

  chomp($s);
  $s =~ s/^\s*//;
  $s =~ s/\s*$//;

  return $s;
}

sub fit_line {

  my $ra = shift;

  my @new = ();
  my $part;
  my $len = 66;

  foreach my $line (@$ra) {
    if( length($line) <= $len ) {
      push(@new,$line);
    } else {
      ($part,$line) = get_part($line,$len);
      while( $line ) {
        push(@new,$part."\\");
        ($part,$line) = get_part($line,$len-16);
	$part = "\t\t" . $part;
      }
      push(@new,$part);
    }
  }

  return \@new;
}

sub get_part {

  my ($line,$len) = @_;

  if( length($line) <= $len ) {
    return ($line,"");
  } else {
    my @words = split(/\s+/,$line);
    my $part = "";
    my $rest = "";
    my $full = 0;
    foreach my $word (@words) {

      if( $full or ( length($part . $word) > $len )) {
	$rest .= "$word ";
	$full = 1;
      } else {
	$part .= "$word ";
      }
    }
    return ($part,$rest);
  }
}


