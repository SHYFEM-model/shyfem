#!/usr/bin/perl -w
#
# changes .f files to .f90 format
#
# changes comment char to !
# changes tabs to spaces in continuation lines
# inserts & at column 73 and 6 for continuation
#
#------------------------------------------

use strict;

my $debug = 0;

my $old_line;
my $line;
my $conti_lines = 0;
my $in_conti = 0;
$::nlines = 0;

while(<>) {

  chomp;

  $::nlines++;

  s/^c/!/;
  s/^C/!/;
  s/^\*/!/;

  if( /^     \S/ ) {
    s/^     \S/     \&/;
    #$_ = insert_ampersand($_) unless $in_conti;
    $old_line = insert_ampersand($old_line);
    $conti_lines++;
    $in_conti = 1;
  } else {
    if( $in_conti ) {
      $old_line = tabs2spaces($old_line);	#this is last conti
    }
    $in_conti = 0;
  }
  
  print "$old_line\n" if defined $old_line;

  $old_line = $_;
}

print "$old_line\n";

print STDERR "$conti_lines continuation lines found\n" if $debug;

#------------------------------------------

sub insert_ampersand {

  my $line = shift;

  #$line .= " &";

  $line = tabs2spaces($line);

  my $l = length($line);

  if( $l > 72 ) {
    my $c73 = substr($line,72,1);
    if( $c73 ne '&' ) {
      die "line $::nlines too long... cannot insert ampersand: $line\n"
    } else {
      return $line;
    }
  }

  $l = 72 - $l;
  while( $l-- ) {
    $line .= " ";
  }
  $line .= "&";

  return $line;
}

sub tabs2spaces {

  # this comes from the perl cookbook

  my $string = shift;

  while ($string =~ s/\t+/' ' x (length($&) * 8 - length($`) % 8)/e) {
    # spin in empty loop until substitution finally fails
  }

  return $string
}

#------------------------------------------

