#!/usr/bin/perl 
####!/usr/bin/perl -i.bak

$param = shift;
$value = shift;

while(<>) {
  if( /^\s+parameter/ ) {
    if( /(\W$param\s*=\s*)([+-]?\d+)/ ) {
        print STDERR "$1 -> $2\n";
    }
  }
}

