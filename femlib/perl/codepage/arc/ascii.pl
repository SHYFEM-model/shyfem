#!/usr/bin/perl -w

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib '/home/georg/lib/perl';
use ascii;

my $ascii = new ascii;

#$ascii->print_actual();

while(<>) {

  chomp;

  my $text = $ascii->decode_text($_);

  print "    $_    ->    $text\n";
}
