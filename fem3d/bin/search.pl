#!/usr/bin/perl

while(<>) {

  #if( /nlv\s*=/ ) {
  if( /hlv\s*\(\s*\S+\s*\)\s*=/ ) {
    print "$ARGV: $_";
  }
}

