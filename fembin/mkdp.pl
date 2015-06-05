#!/usr/bin/perl -ws
#
# checks include and computes dipendencies for fortran and c files
#
# handles recursive includes
#
# 26.02.2015	ggu	adapted for f90 files
#
#------------------------------------------------------

use strict;

$::moddir = "" unless $::moddir;
$::nomake = 1 if $::nomake;
$::debug = 1 if $::debug;

$::make = 1 unless $::nomake;

my $mkdp = "# DO NOT DELETE THIS LINE -- make depend depends on it.";

my @lines = ();

foreach my $file (@ARGV) {

  print STDERR "$file\n" if $::debug;

  my $rlist = {};
  handle_file($file,$rlist,\@lines);

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

sub handle_file {

  my ($file,$rlist,$rlines) = @_;

  my $hfile;
  my $mfile;
  my $fh;
  my %modules_in_file = ();

  open($fh,"$file") || die "Cannot open file $file\n";

  while( <$fh> ) {
    s/\s*\!.*$//;		# get rid of trailing comments
    if( /^\s+include\s*['"]\s*([\w.]+)\s*['"]\s*$/) {
      $hfile = $1;
    } elsif( /^\s*\#\s*include\s*['"]\s*([\w.]+)\s*['"]\s*$/) {
      $hfile = $1;
    } elsif( /^\s+use\s+(\w+)\s*,\s*only\s*:/) {
      $mfile = "$1.mod";
    } elsif( /^\s+use\s+(\w+)\s*,$/) {
      print STDERR "*** cannot handle more than 1 module per line yet\n";
    } elsif( /^\s+use\s+(\w+)\s*$/) {
      $mfile = "$1.mod";
    } elsif( /^\s+module\s+(\w+)\s*$/) {	#must treat differently
      my $module = $1;
      $modules_in_file{"$module.mod"} = 1;
      my $fileo = $file;
      $fileo =~ s/\.f$/.o/;
      $fileo =~ s/\.f90$/.o/;
      my $line = "$module.mod: $fileo";
      $line = "$::moddir/$line" if $::moddir;
      push(@$rlines,$line);
    }

    if( $hfile ) {
      print STDERR "include found: $hfile\n" if $::debug;
      my $ins = insert_inc($rlist,$hfile);
      if( $ins ) {	#new
        handle_file($hfile,$rlist);
      }
      $hfile = "";
    } elsif( $mfile ) {
      print STDERR "module found: $mfile\n" if $::debug;
      if( $modules_in_file{$mfile} ) {
	print STDERR "*** avoid circular reference for module: $mfile\n";
      } else {
        $mfile = "$::moddir/$mfile" if $::moddir;
        insert_inc($rlist,$mfile);
      }
      $mfile = "";
    }
  }

  close($fh);
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

#-----------------------------------------------------------------------

sub write_inc {

  my ($file,$rlist) = @_;
  
  my @list = keys %$rlist;

  return unless scalar @list;

  $file =~ s/\.[fc]$/.o:/;
  $file =~ s/\.f90$/.o:/;
  foreach my $f (@list) {
    $file .= " $f";
  }

  return $file;
}

sub insert_inc {

  my ($rlist,$hfile) = @_;

  if( $rlist->{$hfile} ) {		# already inserted
    return 0;
  } else {
    $rlist->{$hfile}++;
    return 1;
  }
}

#-----------------------------------------------------------------------

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

#-----------------------------------------------------------------------

