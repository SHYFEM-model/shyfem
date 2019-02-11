#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#make_iso_8859_1();
make_windows_1257();

#-------------------------------------------------

sub make_windows_1257 {

# =20     U+0020  SPACE
# =21     U+0021  EXCLAMATION MARK

# 0x80    0x20AC  #EURO SIGN
# 0x82    0x201A  #SINGLE LOW-9 QUOTATION MARK

  while(<>) {

    chomp;
    next unless /^0x/;

    my ($hex,$uni,@rest) = split;

    $hex =~ s/^0x/=/;
    $uni =~ s/^0x/U\+/;
    my $rest = join(" ",@rest);
    $rest =~ s/^\#//;
    $uni = "U+0000" if $uni eq "0";

    print "$hex\t$uni\t$rest\n";
  }

}

sub make_iso_8859_1 {

%entity = ();
@entity = ();

my $what;

while(<>) {

  chomp;
  my @f = split;
  my @new = ();

  push(@new,shift(@f));
  push(@new,shift(@f));
  push(@new,shift(@f));

  $what = shift(@f);

  if( $what =~ /^\&/ or $what =~ /^<em>\&/ ) {
    push(@new,$what);
    push(@new,shift(@f));
  } else {
    push(@new,"--");
    push(@new,$what);
  }

  push(@new,shift(@f));
  push(@new,shift(@f));
  
  $what = shift(@f);

  if( $what ) {
    push(@new,$what);
  } else {
    push(@new,"--");
  }

  my $line = join("  ",@new);
  #print "$line\n";

  insert_entity(\%entity,\@entity,@new[0,1,2,3]);
  insert_entity(\%entity,\@entity,@new[4,5,6,7]);
}

my $i = 0;
foreach my $val (@entity) {
  my $hex = sprintf("%02x",$i);
  print "$i  ($hex)   $val     $entity{$hex}\n";
  $i++;
}

my $i = 0;
foreach my $val (@entity) {
  $aux = quotemeta($val);
  print " \"$aux\"";
  $i++;
  if( $i == 256 ) {
    print "\n";
  } elsif( $i%5 == 0 ) {
    print " ,\n";
  } else {
    print " , ";
  }
}
print "\n";

}

sub insert_entity {

  my $rhash = shift;
  my $rarray = shift;

  my $dec = @_[1];
  my $key = @_[2];
  my $val = @_[3];

  if( $val eq "--" ) {
    $val = @_[0];
  }
  if( $dec < 32 or $dec == 127 ) {
    $val = "--";
  } elsif( $dec == 32 ) {
    $val = " ";
  }
  if( $dec < 16 ) {
    $key = "0" . $key;
  }

  $val =~ s/<em>//g;
  $val =~ s/<\/em>//g;
  $val =~ s/^\&amp;/\&/g;

  $rhash->{$key} = $val;
  push(@$rarray,$val);
}


