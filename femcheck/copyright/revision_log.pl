#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# header
#   header0		empty if no copy, everything in headermain
#   copy		empty lines before and after not counted
#   headermain		everything else
# headermain
#   header1		everything if no revision log
#   rev			start with "revision log :", not excluding this line
#   header2		start with no date line after rev, emtpy if no rev
# body

use lib ("$ENV{SHYFEMDIR}/femcheck/copyright"
		,"$ENV{HOME}/shyfem/femcheck/copyright");

use strict;
use revision_log;

$::file = $ARGV[0];

$::warn = 0;			#warn for revision out of revision log section
$::obsolete = 0;		#check for obsolete date in text

$::extract = 0 unless $::extract;	#extract revlog to revlog.tmp
$::check = 0 unless $::check;		#checks files
$::gitrev = 0 unless $::gitrev;		#uses gitrev for revision log

$::manual = 0;

$::names = ();
make_names();

my @header = ();
my $header = \@header;
my @body = ();
my $body = \@body;

my $in_header = 1;

#print STDERR "reading file $::file\n";

while(<>) {

  chomp;

  if( /^\s*$/ ) {			# empty line
  } elsif( /^[!cC\*]/ ) {		# comment in first col
  } elsif( /^\s*!/ ) {			# comment with leading blanks
  } else {				# first line of code
    $in_header = 0;
  }

  if( $in_header ) {
    push(@$header,$_);
  } else {
    push(@$body,$_);
  }
}

#my $n = @$header;
#print "size of header = $n\n";

my ($header0,$copy,$headermain) = extract_copy($header);
my ($header1,$rev,$header2) = extract_rev($headermain);

my $devs = count_developers($rev);

if( $::gitrev ) {
  if( @$rev == 0 ) {
    $header2 = copy_divisor($header1);
  }
  $rev = integrate_revlog($rev);
}

write_file("$::file.new",$header0,$copy,$header1,$rev,$header2,$body);
#print_file($header0,$copy,$header1,$rev,$header2);
#stats_file($::file,$rev,$devs);

my $has_revision_log = @$rev + $::manual;

exit $has_revision_log;

#--------------------------------------------------------------

sub extract_copy
{
  my $header = shift;

  my $in_copy = 0;
  my @copy = ();
  my @header0 = ();
  my @headermain = ();

  foreach (@$header) {
    #print "$in_copy: $_\n";
    if( /^!------------------------------/ ) {
      $in_copy++;
      if( $in_copy < 3 ) {
        push(@copy,$_);
        next;
      }
    }
    if( $in_copy == 0 ) {
      push(@header0,$_);
    } elsif( $in_copy == 1 ) {
      $::copyright++ if( /^..\s*Copyright/ );
      $::shyfem++ if( /^..\s*This file is part of SHYFEM./ );
      $::manual++ if( /^..\s*This file is part of SHYFEM.\s+\(m\)/ );
      push(@copy,$_);
    } else {
      push(@headermain,$_);
    }
  }

  if( $::copyright ) {
    #if( $::copyright > 1 ) {
    #  print STDERR "*** more than one copyright in file $::file\n";
    #}
    if( $in_copy < 2 ) {
      print STDERR "*** copyright notice not finished in file $::file\n";
    }
    if( $::shyfem == 0 ) {
      print STDERR "*** missing shyfem line in file $::file\n";
    }
  } else {
    print STDERR "*** no copyright in file $::file\n";
    @headermain = @$header;
    @copy = ();
    @header0 = ();
  }

  return (\@header0,\@copy,\@headermain);
}

sub extract_rev
{
  my $headermain = shift;

  my $revs = 0;
  my $in_rev = 0;
  my @rev = ();
  my @header1 = ();
  my @header2 = ();

  foreach (@$headermain) {
    if( /^..\s*revision log :/ ) {
      $in_rev++;
    } elsif( /^[cC!\*]\s*$/ or /^\s*$/ ) {
      if( $in_rev == 1 and $revs > 1 ) {
        $in_rev++;				#account for empty line
      }
    }
    if( $in_rev == 0 ) {
      push(@header1,$_);
    } elsif( $in_rev == 1 ) {
      $revs++;
      push(@rev,$_);
    } else {
      push(@header2,$_);
    }
  }

  if( $in_rev > 2 ) {
    print "in_rev: $in_rev\n";
    print STDERR "*** more than one revision log in file $::file\n";
  }

  return (\@header1,\@rev,\@header2);
}

sub write_file
{
  my $file = shift;

  open(FILE,">$file");

  foreach my $block (@_) {
    foreach my $line (@$block) {
      print FILE "$line\n";
    }
  }

  close(FILE);
}

sub print_file
{
  my $n = 0;

  foreach my $block (@_) {
    $n++;
    print "======================================================\n";
    print "$n\n";
    print "======================================================\n";
    foreach my $line (@$block) {
      print "$line\n";
    }
  }
}

sub stats_file
{
  my ($file,$rev,$dev) =@_;

  my $revs = @$rev;
  $revs -= 2 if $revs > 0;

  my $error = 0;
  $error++ if $::copyright != 1;
  $error++ if $::shyfem != 1;
  $error++ if $::manual != 0;
  $error++ if $revs == 0;

  return unless $error;

  print" $::copyright $::shyfem $::manual $revs $devs   $file\n";
}

sub count_developers
{
  my $rev = shift;

  my @rev = @$rev;
  my $n = @rev;
  return 0 unless $n;

  shift(@rev); shift(@rev);

  foreach (@rev) {
    $n = check_new_revision();
    unless( $n ) {
      print STDERR "*** $::file: cannot parse revlog: $_\n";
    }
  }

  my @keys = keys %::devcount;
  $n = @keys;

  return $n;
}

sub integrate_revlog
{
  my $rev = shift;

  my $nrev = @$rev;
  die "not ready merging, only integrating\n" if $nrev;

  my @gitrev=`git-file -revlog $::file`;
  #my $ngit = @gitrev;
  #print "checking file $::file\n";
  #print "from git $ngit total lines read...\n";

  my @new = ();
  my $in_revision_log = 0;
  foreach (@gitrev) {
    $in_revision_log = 1 if /revision log/;
    next unless $in_revision_log;
    chomp;
    push(@new,$_);
  }
  #pop(@new) if $new[-1] =~ /^\s*$/;

  my $nnew = @new;
  print "from git $nnew lines read...\n";

  return \@new;
}

sub copy_divisor
{
  my $header = shift;

  my @aux = reverse(@$header);
  my @new = ();

  foreach (@aux) {
    push(@new,$_);
    last if /^[cC!]\s*\*\*\*/;
    last if /^[cC!]\s*---/;
    last if /^[cC!]\s*===/;
  }

  @aux = reverse(@new);

  return \@aux;
}

