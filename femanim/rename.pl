#!/usr/bin/perl

while( $file = shift ) {
  @f = split(/\./,$file);

  $n = $f[1];

  if( $n < 10 ) {
    $num = "00" . $n;
  } elsif( $n < 100 ) {
    $num = "0" . $n;
  } else {
    $num = $n;
  }

  $f[1] = $num;

  $newfile = join(".",@f);

  print "$file -> $newfile\n";

  rename $file , $newfile;
}
