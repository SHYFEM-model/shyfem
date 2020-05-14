#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

bin/newstruct.pl -struct=ht *.f > g1.txt	# create structure

cat g1.txt | tr -d " \t" > g2.txt		# delete spaces

sort g2.txt | uniq  -c | sort -n > g3.txt	# find occurence of calls

cat g3.txt

