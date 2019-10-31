#!/usr/bin/perl

while(<>) {

  if( /define LOG_INFO/ ) {
    print; print "    \\\n";
    <>;
    next;
  }

  print;
}

