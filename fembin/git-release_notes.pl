#!/usr/bin/perl -s

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

$debug = 0;
$title = "" unless $title;

$tag = shift;

unless( $tag ) {
  die "no tag given... Usage: git-release_notes.pl tag RELEASE_NOTES\n";
}
unless( $ARGV[0] ) {
  die "no file given... Usage: git-release_notes.pl tag RELEASE_NOTES\n";
}

print STDERR "=== using tag $tag with title $title\n" if $debug;

$text = "\n";
$title_text = "";
$in_note = 0;

while(<>) {

  if( /\s+$tag\s*(.*)\s*$/ ) {
    die "more than one entry with this tag: $tag\n" if $in_note;
    $in_note = 1;
    $title_text = $1;
    print STDERR "=== extracting title: $title $title_text\n" if $debug;
  } elsif( /VERS_/ ) {
    last if $in_note;
  }

  if( $in_note ) {
    $text .= $_;
  }

}

if( $title ) {
  print "$title_text";
} else {
  print "$text";
}

