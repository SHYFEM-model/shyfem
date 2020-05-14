#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# automerges merge conflicts
#
#-----------------------------------------------------------

my $file = $ARGV[0];

$debug = 0;

$start="<<<<<<< HEAD";
  $end=">>>>>>>";
  $mid="=======";

  $sep="============================================================";

print STDERR "handling file $file\n" if $debug;

while(<>) {

  if( /^$start$/ ) {
    $block++;
    $in_first = 1;
    next;
  } elsif( /^$mid$/ ) {
    if( not $in_first ) {
      die "*** error reading mid but not in first: $_";
    }
    $in_first = 0;
    $in_second = 1;
    next;
  } elsif( /^$end / ) {
    if( not $in_second ) {
      die "*** error reading end but not in second: $_";
    }
    $in_second = 0;
    handle_blocks();
    next;
  }

  if( $in_first ) {
    push(@first,$_);
  } elsif( $in_second ) {
    push(@second,$_);
  } else {
    print;
  }
}

print STDERR "finished handling file $file with $block blocks(s)\n" if $debug;

#-----------------------------------------------------------

sub handle_blocks
{
  my $ftime = get_time(@first);
  my $stime = get_time(@second);

  #print STDERR "first: $ftime   second: $stime\n";

  if( $ftime < $stime ) {
    my @aux = @second; @second = @first; @first = @aux;
  }

  print_block(@first);
  print_sep();
  print_block(@second);

  @first = (); @second = ();
}

sub print_block
{
  foreach (@_) {
    print;
  }
}

sub print_sep
{
  if( $file eq "VERSION" ) {
    die "*** can handle only one block\n" if( $block > 1 );
    print "\n";
  } elsif( $file eq "COMMIT" ) {
    die "*** can handle only one block\n" if( $block > 1 );
    print "\n$sep\n\n";
  } elsif( $file eq "fem3d/subver.f" ) {
    @second = ();	#we only need the first occurrence, no seperator
  } else {
    die "*** cannot handle file: $file\n";
  }
}

sub get_time
{
  my $date = shift;

  chomp($date);

  if( $file eq "VERSION" ) {
    $date =~ s/^.*commit_//;
  } elsif( $file eq "COMMIT" ) {
    ;
  } elsif( $file eq "fem3d/subver.f" ) {
    $date =~ s/^.*(\d{4}-\d{2}-\d{2}).\s*$/$1/;
  } else {
    die "*** cannot handle file: $file\n";
  }

  print STDERR "date found: $date\n" if $debug;
  my $secs = `date -d "$date" +%s`;

  return $secs;
}

#-----------------------------------------------------------

