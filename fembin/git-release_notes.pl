#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

$tag = shift;

unless( $tag ) {
  die "no tag given... Usage: git-release_notes.pl tag RELEASE_NOTES\n";
}
unless( $ARGV[0] ) {
  die "no file given... Usage: git-release_notes.pl tag RELEASE_NOTES\n";
}

#print STDERR "using tag $tag\n";

$text = "\n";
$in_note = 0;

while(<>) {

  if( /\s+$tag$/ ) {
    die "more than one entry with this tag: $tag\n" if $in_note;
    $in_note = 1;
  } elsif( /VERS_/ ) {
    last if $in_note;
  }

  if( $in_note ) {
    $text .= $_;
  }

}

print "$text";

