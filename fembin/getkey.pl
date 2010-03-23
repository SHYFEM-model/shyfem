#!/usr/bin/perl -s

if( $h or $ $help or not $ARGV[0] ) {
  die "Usage: getkey.pl keyword file\n";
}

$keyword = shift;
$key = quotemeta($keyword);

#print STDERR "keyword: $keyword - $key\n";

while(<>) {

  if( s/^\s*$key\s*:?\s*// ) {
    print;
  }

}

