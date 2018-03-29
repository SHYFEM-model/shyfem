#!/usr/bin/perl
#
# splits xrange in lower and upper part
#
# year is needed, other parts can be obmitted
# format: YYYY-MM-DD::hh:mm:ss
#
#---------------------------------------------

$xrange = $ARGV[0];

#print STDERR "xrange (perl): $xrange\n";

$xhigh = "";
$xlow = "";

if( $xrange =~ /^\s*$/ ) {
  ;
} elsif( $xrange =~ /^:(.*)$/ ) {		# only high
  $xhigh = regular($1);
} elsif( $xrange =~ /^(.*):$/ ) {		# only low
  $xlow = regular($1);
} elsif( $xrange =~ /^(.*):(.*)::(.*)$/ ) {	# :: is from high
  $xlow = regular($1);
  $xhigh = regular("$2::$3");
} elsif( $xrange =~ /^(.*):(.*)-(.*)$/ ) {	# - is in date of high
  $xlow = regular($1);
  $xhigh = regular("$2-$3");
} else {
  die "cannot parse xrange: $xrange\n";
}

print "$xlow=$xhigh\n";

#---------------------------------------------

sub regular {

  my $string = shift;

  my ($dd,$tt);
  my ($y,$m,$d);
  my ($h,$M,$s);
  my ($out);

  if( $string =~ "::" ) {
    ($dd,$tt) = split(/::/,$string);
  } elsif( $string =~ /:/ ) {
    $tt = $string;
  } elsif( $string =~ /-/ ) {
    $dd = $string;
  } elsif( $string =~ /^(\d\d\d\d)$/ ) {
    $dd = $string;
  } else {
    die "cannot regularize date (something is missing): $string\n";
  }

  if( not $dd ) {
    die "missing date in string: $string\n";
  } elsif( $dd =~ /^(\d\d\d\d)$/ ) {
    $y = $1;
  } elsif( $dd =~ /^(\d\d\d\d)-(\d\d)$/ ) {
    $y = $1;
    $m = $2;
  } elsif( $dd =~ /^(\d\d\d\d)-(\d\d)-(\d\d)$/ ) {
    $y = $1;
    $m = $2;
    $d = $3;
  } else {
    die "cannot parse date: $dd ($string)\n";
  }

  $m = "01" unless $m;
  $d = "01" unless $d;

  if( not $tt ) {
    ;
  } elsif( $tt =~ /^(\d\d)$/ ) {
    $h = $1;
  } elsif( $tt =~ /^(\d\d):(\d\d)$/ ) {
    $h = $1;
    $M = $2;
  } elsif( $tt =~ /^(\d\d):(\d\d):(\d\d)$/ ) {
    $h = $1;
    $M = $2;
    $s = $3;
  } else {
    die "cannot parse time: $tt ($string)\n";
  }

  $h = "00" unless $h;
  $M = "00" unless $M;
  $s = "00" unless $s;

  #print STDERR "regular: $string $dd $tt\n";
  #print STDERR "regular: $y $m $d :: $h $M $s\n";

  return "$y-$m-$d\:\:$h:$M:$s";
}

#---------------------------------------------

