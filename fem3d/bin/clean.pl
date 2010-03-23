#!/usr/bin/env perl

%trans = (
		"10"	=>	"\n",		#(0A) -> \n
		"13"	=>	"",		#(0D) -> \r
		"146"	=>	"'",		#(92) -> '
		"150"	=>	"-",		#(96) -> -  (?)
		"133"	=>	" ",		#(85) -> 'space'
	 );

while(<>) {

  $n = length($_);

  for( $i=0 ; $i < $n ; $i++ ) {
    $c = substr($_,$i,1);
    $nc = ord($c);
    if( $nc < 32 || $nc > 127 ) {
      $count[$nc]++;
      $newc = &subst($nc);
      #print "$nc\n";
      $c = $newc;
    }
    print $c;
  }

}

$n = @count;
$count[10] = 0;
for( $i=0 ; $i < $n ; $i++ ) {
  if( $count[$i] ) {
    $hnc = &tohex($i);
    print STDERR "$i ($hnc):  $count[$i]\n";
  }
}

sub tohex {

  my $nc = shift;
  my $hex = sprintf("%02.2X",$nc);
  return $hex;
}

sub subst {

  my $nc = shift;

  if( exists $trans{$nc} ) {
    return $trans{$nc};
  } else {
    $hnc = &tohex($nc);
    print STDERR "Cannot convert: $nc  ($hnc)\n";
    return chr($nc);
  }
}
  
