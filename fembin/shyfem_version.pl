#!/usr/bin/perl -s
#
# parses VERSION file
#
# Usage: shyfem_version.pl [option] VERSION
#
#    options:
#		-version
#		-date
#		-tag
#		-extra
#		-noextra
#		-tag_extra
#
#-----------------------------------------------------

$_ = <>;
chomp;

@f = split;
$n = @f;

die "cannot parse first line of VERSION file:\n$_\n" unless /^version/;

if( $noextra ) {
  s/\s+\S+\s*$// if $n >= 5;
} elsif( $version ) {
  $_ = $f[1];
} elsif( $date ) {
  $_ = $f[2];
} elsif( $tag ) {
  $_ = $f[3];
} elsif( $extra ) {
  $_ = $f[4];
} elsif( $tag_extra ) {
  $_ = $f[3];
  $_ .= "_$f[4]" if $n >= 5;
}

print "$_\n";

