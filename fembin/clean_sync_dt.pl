#!/usr/bin/perl
#
# deletes lines refering to syncronization time step
#
#------------------------------------------------------

while(<>) {

  @f = split;
  next if $f[3];		# 1 indicates syncronization time step
  print;
}

#------------------------------------------------------

