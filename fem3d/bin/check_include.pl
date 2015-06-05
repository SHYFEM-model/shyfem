#!/usr/bin/perl
#
#----------------------------

my %inc = ();
my %files = ();

while( <> ) {

  chomp;
  s/\s+//g;
  if( /include[\'\"](\w+)\b/ ) {
    $inc{$1}++;
    $files{$1} .= "$ARGV ";
  }
}

@keys = sort { $inc{$a} <=> $inc{$b} }  keys %inc;

foreach my $inc (@keys) {
  my $count = $inc{$inc};
  my $fc = file_count($files{$inc});
  print "$count   ($fc)   $inc\n";
}

#-------------------------------------

sub file_count {

  my $string = shift;

  my @f = split(/\s+/,$string);

  my %f = ();
  foreach my $file (@f) {
    $f{$file}++;
  }

  my @ff = keys %f;
  my $n = scalar @ff;

  return $n;
}

#-------------------------------------

