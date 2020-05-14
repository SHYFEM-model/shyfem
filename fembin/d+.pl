#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# shows how much diffs

$bin = 10;	#bin size

while( <> ) {

  if( /^d.* (\d+)\n$/ ) {
    $del += $1;
  } elsif( /^a.* (\d+)\n$/ ) {
    $add += $1;
    &skip_lines($1);
  } else {
    die "Internal error : $_\n";
  }
}

#print "$add/$del ";

$total = $add + $del;
$total = &format_number($total,5);

while( $add > 0 ) {
  $a .= "+";
  $add -= $bin;
}

while( $del > 0 ) {
  $d .= "-";
  $del -= $bin;
}

print "$total $a$d\n";

#######################

sub skip_lines {

  my $n = $_[0];

  while( $n ) {
    $_ = <>;
    $n--;
  }
}

sub format_number {

  my $num = shift;
  my $dig = shift;

  my $new = "";
  my $i;
  my $l = length($num);

  for($i=0;$i<$dig-$l;$i++) {
    $new .= " ";
  }

  return $new . $num;
}
