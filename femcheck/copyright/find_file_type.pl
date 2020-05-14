#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femcheck/copyright"
		,"$ENV{HOME}/shyfem/femcheck/copyright");

use strict;
use revision_log;

my $file = $ARGV[0];
my $type = find_file_type($file);

print "$type\n";

#------------------------------------------------------------------------

