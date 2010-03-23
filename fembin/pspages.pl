#!/usr/bin/perl

while(<>) {

  if( /^%%Pages:\s+(\d+)/ ) {
    $pages = $1;
    die "$pages\n";
  }
}

print "0\n";

