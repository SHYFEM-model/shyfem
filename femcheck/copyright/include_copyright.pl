#!/usr/bin/perl -w
#
#-----------------------------------

use strict;

my $copyright = "/home/georg/shyfem/femcheck/copyright/copyright_notice.txt";

my $lang = shift;

if( $lang eq "f" ) {
  ;
} elsif( $lang eq "c" ) {
  ;
} else {
  print "Usage: include_copyright.pl language file\n";
  die "unknown language: $lang ... please choose between f or c\n";
}

my $copy = read_copyright();

print "@$copy";

$copy = adjust_copyright($copy,$lang);

print "@$copy";








#-----------------------------------

sub adjust_copyright
{
  my ($copy,$lang) = @_;

  my ($char,$bline,$eline);

  if( $lang eq "f" ) {
    $char = "!";
    $bline = "\
!--------------------------------------------------------------------------\n";
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
    push(@new,"$char$_");
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

