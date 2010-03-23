#!/usr/bin/perl -w

use lib '/home/georg/lib/perl';
use ascii;

my $ascii = new ascii;

#$ascii->print_actual();

while(<>) {

  chomp;

  my $text = $ascii->decode_text($_);

  print "    $_    ->    $text\n";
}
