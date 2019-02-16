#!/usr/bin/perl -ws
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# checks include and computes dipendencies for fortran and c files
#
# handles recursive includes
#
# 26.02.2015	ggu	adapted for f90 files
# 04.01.2016	ggu	sort targets and dependencies alphabetically
#
#------------------------------------------------------

use strict;

#------------------------------------------------------

%::ignore_modules = (
			'mpi.mod'	=> 	1,
		    );

#------------------------------------------------------

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
  print "$line\n" if ($line and $::debug);
  push(@lines,$line) if $line;
}

if( not $::make ) {
  print STDERR "Makefile not changed...\n";
  exit 0;
}

my $mfile = change_makefile(\@lines);

if( is_different("$mfile","$mfile.new") ) {
  print STDERR "Makefile changed... substituting\n";
  rename("$mfile","$mfile.bak");
  rename("$mfile.new","$mfile");
} else {
  print STDERR "Makefile did not change... not substituting\n";
  rename("$mfile.new","$mfile.bak");
}

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
    if( /^\s+include\s*['"]\s*([\w.]+)\s*['"]\s*$/i) {
      $hfile = $1;
    } elsif( /^\s*\#\s*include\s*['"]\s*([\w.]+)\s*['"]\s*$/i) {
      $hfile = $1;
    } elsif( /^\s*use\s+(\w+)\s*,\s*only\s*:/i) {
      my $module = lc($1);
      $mfile = "$module.mod";
    } elsif( /^\s*use\s+(\w+)\s*,$/i) {
      print STDERR "*** cannot handle more than 1 module per line yet\n";
    } elsif( /^\s*use\s+(\w+)\s*$/i) {
      my $module = lc($1);
      $mfile = "$module.mod";
    } elsif( /^\s*module\s+(\w+)\s*$/i) {	#must treat differently
      my $module = lc($1);
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
      if( $::ignore_modules{$mfile} ) {
	if( $::debug ) {
	  print STDERR "*** ignoring module: $mfile\n";
	}
      } elsif( $modules_in_file{$mfile} ) {
	if( $::debug ) {
	  print STDERR "*** avoid circular reference for module: $mfile\n";
	}
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
        print STDERR "make line found...\n" if $::debug;
	$mkdp_found = 1;
        last;
    }
    print NEW;
  }

  print NEW "\n" unless $mkdp_found;		# just for beauty
  print NEW "$mkdp\n";
  print NEW "\n";

  #$ra = fit_line($ra);
  my @list = sort(@$ra);
  $ra = fit_line(\@list);

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
  
  #my @list = keys %$rlist;
  my @list = sort keys %$rlist;

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

sub is_different {

  my ($f1,$f2) = @_;

  system("cmp --quiet $f1 $f2");
  my $status = $?;
  #print STDERR "exit code is $status\n";

  return $status;
}

#-----------------------------------------------------------------------

