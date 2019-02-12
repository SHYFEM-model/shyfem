#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

@text = <>;
$file = $ARGV;
$name = $file;
$name =~ s/\.tmp//;
$name = "\n" . "_" . $name . "_data = [";

$first = shift(@text);
$first = "$name" . $first;
unshift(@text,$first);

$last = pop(@text);
$last =~ s/\],/\]\]\n/;
push(@text,$last);

print @text;

