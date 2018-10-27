#!/usr/bin/perl -w
#
#-----------------------------------

use strict;

my $copyright = "/home/georg/shyfem/femcheck/copyright/copyright_notice.txt";

my $lang = shift;
$lang = "" unless $lang;

if( $lang eq "f" ) {
  ;
} elsif( $lang eq "c" ) {
  ;
} else {
  print "Usage: include_copyright.pl language file\n";
  die "unknown language: $lang ... please choose between f or c\n";
}

my $copy = read_copyright();
$copy = adjust_copyright($copy,$lang);
print_copyright($copy);

print_file();

#-----------------------------------

sub adjust_copyright
{
  my ($copy,$lang) = @_;

  my ($char,$bline,$eline,$line);

  if( $lang eq "f" ) {
    $char = "!";
    $bline = 
"!--------------------------------------------------------------------------\n";
    $eline = $bline;
  } else {
    $char = "*";
    $bline = "/************************************************************\n";
    $eline = "************************************************************/\n";
  }

  my @new = ();

  push(@new,"\n");
  push(@new,$bline);

  foreach (@$copy) {
    $line = $char . $_;
    push(@new,$line);
  }

  push(@new,$eline);
  push(@new,"\n");

  return \@new;
}

sub read_copyright
{
  open(COPY,"<$copyright") || die "cannot read file: $copyright\n";
  my @copy = <COPY>;
  close(COPY);

  return \@copy;
}

sub print_copyright
{
  my $a = shift;

  foreach (@$a) {
    print "$_";
  }
}

sub print_file
{
  my $debug = 0;	# set to > 0 for debug
  my $header = 1;	# we are still in header

  while(<>) {
    if( $header ) {
      next if /^[cC!]\s*$/;
      next if /^\s*$/;
      next if /^[cC!]\s+\$Id:/;
    }
    last if( --$debug == 0 );	# only executed if debug was > 0
    print;
    $header = 0;
  }
}

