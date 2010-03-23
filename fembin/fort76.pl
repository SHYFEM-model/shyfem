#!/usr/bin/perl
#
# extracts single points from fort.76 file

@f = &readline;
$nodes = $f[0];

# skip other header

for($i=0;$i<$nodes;$i++) {
  <>;
}

# open files

for($i=1;$i<=$nodes;$i++) {
  $file = "z.$i"; $filehandle = "z$i";
  open($filehandle,">$file");
  $file = "u.$i"; $filehandle = "u$i";
  open($filehandle,">$file");
  $file = "v.$i"; $filehandle = "v$i";
  open($filehandle,">$file");
  $file = "s.$i"; $filehandle = "s$i";
  open($filehandle,">$file");
}

# loop through file

while ( @f = &readline ) {

  $time = $f[0];

  for($i=1;$i<=$nodes;$i++) {

    ($n,$u,$v,$z,$dummy) = &readline;

    $filehandle = "z$i";
    print $filehandle "$time $z\n";
    $filehandle = "u$i";
    print $filehandle "$time $u\n";
    $filehandle = "v$i";
    print $filehandle "$time $v\n";
    $filehandle = "s$i";
    $s = sqrt( $u*$u + $v*$v );
    print $filehandle "$time $s\n";
  }

}

# close files

for($i=1;$i<=$nodes;$i++) {
    $filehandle = "z$i";
    close( $filehandle );
    $filehandle = "u$i";
    close( $filehandle );
    $filehandle = "v$i";
    close( $filehandle );
    $filehandle = "s$i";
    close( $filehandle );
}

###################################################################

sub readline {

  $_ = <>;

  s/^\s+//;
  s/\s*\n$//;

  return split;

}
