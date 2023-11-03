#!/usr/bin/perl
#
# changes .f files to .f90 format
#
#------------------------------------------

while(<>) {

  chomp;

  s/^c/!/;
  s/^C/!/;
  s/^\*/!/;

  if( /^     \S/ ) {
    s/^     \S/     \&/;
    $line .= " &";
  }
  
  print "$line\n" if $has_line;

  $has_line = 1;
  $line = $_;
}

print "$line\n";

#------------------------------------------

