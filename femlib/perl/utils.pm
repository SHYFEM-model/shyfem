#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# general utilities
#
# useage:
#
# require "utils.pl";		# ts_util.pl must be in same dir
#
#------------------------------------------------------------------------

use strict;

#------------------------------------------------------------------------

sub example_sort_hash_table
{

    my %planets;

    foreach my $distance (sort {$a <=> $b} values %planets) {
        say $distance;
    }

    foreach my $name (sort { $planets{$a} <=> $planets{$b} } keys %planets) {
        printf "%-8s %s\n", $name, $planets{$name};
    }

    foreach my $name (sort { $planets{$a} <=> $planets{$b} or $a cmp $b } 
							keys %planets) {
        printf "%-8s %s\n", $name, $planets{$name};
    }
}

#------------------------------------------------------------------------
1;
#------------------------------------------------------------------------

