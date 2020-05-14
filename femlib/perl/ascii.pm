#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# ascii utilities
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
# use ascii;
# 
# my $ascii = new ascii;
#
##############################################################

use strict;

package ascii;

##############################################################

sub new
{
    my $self;

    $self =	{
	    		 fall_back	=>	0
	    		,binary_encode	=>	0
	    		,act_encode	=>	"q"
	    		,codepage_dir	=>	"/home/georg/lib/perl/codepage"
	    		,actual		=>	undef
	    		,standard	=>	undef
			,iso_8859_1	=>	{}
		};

    bless $self;
    $self->init_arrays();
    $self->init_modules();
    $self->{actual} = $self->{iso_8859_1};
    $self->{standard} = $self->{iso_8859_1};
    return $self;
}

###############################################################################

sub set_charset
{
	my ($self,$charset) = @_;

	$charset = lc($charset);

	my $return_code = 1;

	if( $charset eq "iso-8859-1" ) {
	  $self->{actual} = $self->{iso_8859_1};
	} elsif( $charset eq "us-ascii" ) {
	  $self->{actual} = $self->{iso_8859_1};
	} elsif( $charset eq "utf-8" ) {
	  $self->{actual} = $self->{iso_8859_1};
	} elsif( $charset eq "unknown" ) {
	  $self->{actual} = $self->{iso_8859_1};
	} else {
	  $return_code = $self->try_charset($charset);
	}

	return $return_code;
}

sub reset_charset
{
	my ($self) = @_;

	$self->{actual} = $self->{standard};
}

sub try_charset
{
	my ($self,$charset) = @_;

	my $return_code = 1;
	my $enc = $charset;
	$enc =~ s/-/_/g;

	if( $self->{$enc} ) {
	  #print STDERR "setting charset $charset ($enc)\n";
	  $self->{actual} = $self->{$enc};
	} elsif( $self->{fall_back} ) {
	  $self->{actual} = $self->{standard};
	} else {
	  $return_code = 0;
	}

	return $return_code;
}

sub set_encoding
{
	my ($self,$encoding) = @_;

	my $return_code = 1;

	if( $encoding eq "b" ) {
	  if( $self->{binary_encode} == 0 ) {
	    $return_code = 0;
	  }
	} elsif( $encoding eq "q" ) {
	  ;	# ok, use quoted printable
	} else {
	  print STDERR "*** unknown encoding = $encoding\n";
	  $return_code = 0;
	}

	if( $return_code == 0 ) {
	  $self->{act_encode} = $encoding;
	}

	return $return_code;
}

sub reset_encoding
{
	my ($self) = @_;

	$self->{act_encode} = "q";
}

#-----------------------------------------------------------------

sub   set_fall_back { $_[0]->{fall_back} = 1; }
sub reset_fall_back { $_[0]->{fall_back} = 0; }

#-----------------------------------------------------------------

sub is_encoded
{
	my ($self,$text) = @_;

	if( $text =~ /=\?/ ) {
	  return 1;
	} else {
	  return 0;
	}
}

sub decode_text
{
	my ($self,$text) = @_;
	my $rhash = $self->{actual};

	my $error_text = "???";

	my ($pre,$post,$charset,$encoding,$chars);

	while( $text =~ /=\?([\w-]+)\?(\w)\?([^?]*)\?=/ ) {
	  $pre = $`;
	  $post = $';
	  $charset = lc($1);
	  $encoding = lc($2);
	  $chars = $3;

	  #print STDERR "|$text|$pre|$post|$charset|$encoding|$chars|\n";

	  if( not $self->set_encoding($encoding) ) {
	    print STDERR "** cannot decode with encoding = $encoding\n";
	    return $error_text;
	  }
	  if( not $self->set_charset($charset) ) {
	    print STDERR "** cannot decode with charset = $charset\n";
	    return $error_text;
	  }

	  if( $encoding eq "b" ) {
	    $chars = $self->convert_from_binary($chars);
	    if( $chars =~ /^=\w\w=\w\w/ ) {
	      print STDERR "** binary string still encoded...\n";
	      $chars = $self->convert_from_quoted($chars);
	    }
	  } else {
	    $chars = $self->convert_from_quoted($chars);
	  }

	  $text = $pre . $chars . $post;
	}

	$self->reset_charset();

	return $text
}

sub convert_from_quoted
{
	my ($self,$chars) = @_;

	while( $chars =~ /(=\w{2})/ ) {
	  my $char = $self->get_char($1);
	  $chars =~ s/=\w{2}/$char/;
	}

	return $chars;
}

sub convert_from_binary
{
	my ($self,$chars) = @_;

	my $decoded = decode_base64($chars);
	my @f = split("",$decoded);

	my $newline;
	foreach my $char (@f) {
	  my $dec = ord($char);
	  my $hex = sprintf("%02x",$dec);
	  $newline .= "=$hex";
	}

	return $newline;
}

sub get_char
{
	my ($self,$encode) = @_;
	my $rhash = $self->{actual};

	$encode =~ s/^=//;
	return $rhash->{lc($encode)};
}

sub print_actual
{
	my ($self) = @_;
	my $rhash = $self->{actual};

	my $i = 0;
	for($i=0;$i<256;$i++) {
	  my $hex = sprintf("%02x",$i);
	  my $val = $rhash->{$hex};
	  print "$i  ($hex)   $val\n";
	}
}

sub make_hash
{
	my ($self,$chars) = @_;
	my %hash = ();

	my $i = 0;
	foreach my $char (@$chars) {
	  my $hex = sprintf("%02x",$i);
	  $hash{$hex} = $char;
	  $i++;
	}

	return \%hash;
}

#---------------------------------------------------------------

sub make_codetable
{
	my ($self,$charset) = @_;

	my $rhash = $self->read_codetable($charset);
	return unless $rhash;

	my $base = $self->{iso_8859_1};
	my %newhash = %$base;		#clone base charset

	foreach my $key (keys %$rhash) {
	  #print STDERR "inserting new key: $key\n";
	  $newhash{$key} = $rhash->{$key};
	}

	my $enc = $charset;
	$enc =~ s/-/_/g;
	$self->{$enc} = \%newhash;
	print STDERR "ascii.pm: new code table ready $charset ($enc)\n";
}

sub read_codetable
{
	my ($self,$file) = @_;

	my $dir = $self->{codepage_dir};

	my $filename = "$dir/$file.txt";
	my %hash = ();

	if( not open(CODE,"<$filename") ) {
	  print STDERR "Cannot open open code table $filename\n";
	  return;
	}

	while(<CODE>) {
	  chomp;
	  my ($hex,$uni,$rest) = split;
	  $hex =~ s/^=//;
	  $hex = lc($hex);
	  $uni =~ s/^U\+//;
	  my $aux = lc(substr($uni,2));
	  if( $hex ne $aux ) {
	    $hash{$hex} = "\&\#" . $uni . "\;";
	  }
	}

	close(CODE);

	return \%hash;
}

sub init_modules
{
	my ($self) = @_;

	my $mod = "MIME::Base64";
	#my $mod = "ggu_unknown";
	eval "use $mod";
	if ( $@ ) {
	  print STDERR "*** cannot load module $mod\n";
	  $self->{binary_encode} = 0;
	} else {
	  $self->{binary_encode} = 1;
	}
}

sub init_arrays
{
	my ($self) = @_;

	$self->{iso_8859_1} = $self->make_iso_8859_1();

	$self->make_codetable("iso-8859-2");
	$self->make_codetable("iso-8859-4");
	$self->make_codetable("iso-8859-9");
	$self->make_codetable("windows-1252");
	$self->make_codetable("windows-1257");
	$self->make_codetable("koi8-r");
}

sub make_iso_8859_1
{
	my ($self) = @_;

	my @chars = (
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\ " ,  "\!" ,  "\&quot\;" ,
 "\#" ,  "\$" ,  "\%" ,  "\&amp\;" ,  "\'" ,
 "\(" ,  "\)" ,  "\*" ,  "\+" ,  "\," ,
 "\-" ,  "\." ,  "\/" ,  "0" ,  "1" ,
 "2" ,  "3" ,  "4" ,  "5" ,  "6" ,
 "7" ,  "8" ,  "9" ,  "\:" ,  "\;" ,
 "\&lt\;" ,  "\=" ,  "\&gt\;" ,  "\?" ,  "\@" ,
 "A" ,  "B" ,  "C" ,  "D" ,  "E" ,
 "F" ,  "G" ,  "H" ,  "I" ,  "J" ,
 "K" ,  "L" ,  "M" ,  "N" ,  "O" ,
 "P" ,  "Q" ,  "R" ,  "S" ,  "T" ,
 "U" ,  "V" ,  "W" ,  "X" ,  "Y" ,
 "Z" ,  "\[" ,  "\\" ,  "\]" ,  "\^" ,
 "_" ,  "\`" ,  "a" ,  "b" ,  "c" ,
 "d" ,  "e" ,  "f" ,  "g" ,  "h" ,
 "i" ,  "j" ,  "k" ,  "l" ,  "m" ,
 "n" ,  "o" ,  "p" ,  "q" ,  "r" ,
 "s" ,  "t" ,  "u" ,  "v" ,  "w" ,
 "x" ,  "y" ,  "z" ,  "\{" ,  "\|" ,
 "\}" ,  "\~" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,  "\-\-" ,
 "\&nbsp\;" ,  "\&iexcl\;" ,  "\&cent\;" ,  "\&pound\;" ,  "\&curren\;" ,
 "\&yen\;" ,  "\&brvbar\;" ,  "\&sect\;" ,  "\&uml\;" ,  "\&copy\;" ,
 "\&ordf\;" ,  "\&laquo\;" ,  "\&not\;" ,  "\&shy\;" ,  "\&reg\;" ,
 "\&macr\;" ,  "\&deg\;" ,  "\&plusmn\;" ,  "\&sup2\;" ,  "\&sup3\;" ,
 "\&acute\;" ,  "\&micro\;" ,  "\&para\;" ,  "\&middot\;" ,  "\&cedil\;" ,
 "\&sup1\;" ,  "\&ordm\;" ,  "\&raquo\;" ,  "\&frac14\;" ,  "\&frac12\;" ,
 "\&frac34\;" ,  "\&iquest\;" ,  "\&Agrave\;" ,  "\&Aacute\;" ,  "\&Acirc\;" ,
 "\&Atilde\;" ,  "\&Auml\;" ,  "\&Aring\;" ,  "\&AElig\;" ,  "\&Ccedil\;" ,
 "\&Egrave\;" ,  "\&Eacute\;" ,  "\&Ecirc\;" ,  "\&Euml\;" ,  "\&Igrave\;" ,
 "\&Iacute\;" ,  "\&Icirc\;" ,  "\&Iuml\;" ,  "\&ETH\;" ,  "\&Ntilde\;" ,
 "\&Ograve\;" ,  "\&Oacute\;" ,  "\&Ocirc\;" ,  "\&Otilde\;" ,  "\&Ouml\;" ,
 "\&times\;" ,  "\&Oslash\;" ,  "\&Ugrave\;" ,  "\&Uacute\;" ,  "\&Ucirc\;" ,
 "\&Uuml\;" ,  "\&Yacute\;" ,  "\&THORN\;" ,  "\&szlig\;" ,  "\&agrave\;" ,
 "\&aacute\;" ,  "\&acirc\;" ,  "\&atilde\;" ,  "\&auml\;" ,  "\&aring\;" ,
 "\&aelig\;" ,  "\&ccedil\;" ,  "\&egrave\;" ,  "\&eacute\;" ,  "\&ecirc\;" ,
 "\&euml\;" ,  "\&igrave\;" ,  "\&iacute\;" ,  "\&icirc\;" ,  "\&iuml\;" ,
 "\&eth\;" ,  "\&ntilde\;" ,  "\&ograve\;" ,  "\&oacute\;" ,  "\&ocirc\;" ,
 "\&otilde\;" ,  "\&ouml\;" ,  "\&divide\;" ,  "\&oslash\;" ,  "\&ugrave\;" ,
 "\&uacute\;" ,  "\&ucirc\;" ,  "\&uuml\;" ,  "\&yacute\;" ,  "\&thorn\;" ,
 "\&yuml\;"
	);

	my $hash = $self->make_hash(\@chars);

	return $hash;
}

################################
1;
################################
