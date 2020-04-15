#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# mail utility routines
#
# Usage:
#
# #!/usr/bin/env perl
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use mail;
# 
# $file = @ARGV[0];
# 
# my $mbox = new mail;
# $mbox->open($file);
# $mbox->info;
#
# or more compact
#
# my $mbox = new mail($file);
# $mbox->info;
#
#
##############################################################

use strict;
#require "date.pl";

package message;

sub new
{
    my ($pck,$env,$header,$body) = @_;

    my $self;

    $self =	{
	    		 envelop	=>	$env
			,body		=>	$body
			,headlist	=>	undef
			,headhash	=>	{}
			,headkeys	=>	[]
		};

    bless $self;

    $self->make_hash_header($header);

    return $self;
}

sub get_as_string
{
    my ($self) = @_;

    my @all;

    push(@all,$self->{envelop});
    push(@all,@{$self->{headlist}});
    push(@all,"\n");
    push(@all,@{$self->{body}});

    my $all = join("",@all);

    return $all;
}

sub get_envelop
{
    my ($self) = @_;

    return $self->{envelop};
}

sub get_body
{
    my ($self) = @_;

    return $self->{body};
}

sub get_headers
{
    my ($self,$key) = @_;

    my $hash = $self->{headhash};
    my $value = $hash->{$key};

    #my @values = split(/\s*,\s*/,$value);
    #my @values = split(/\s*[,;]\s*/,$value);
    my @values = $self->split_addresses($value);

    return @values;
}

sub get_header
{
    my ($self,$key) = @_;

    my $hash = $self->{headhash};

    if( $hash->{$key} ) {
      return $hash->{$key};
    } else {
      return "";
    }
}

sub get_header_first		#gets first header from list of headers
{
    my ($self,@keys) = @_;

    my $hash = $self->{headhash};

    foreach my $key (@keys) {
      my $value = $hash->{$key};
      return $value if defined $value;
    }

    return "";
}

sub append
{
    my ($self,$folder) = @_;

    unless( open(FOLDER,">>$folder") ) {
      print STDERR "Cannot open folder $folder\n";
      return;
    }

    select(FOLDER);
    $self->write();
    select(STDOUT);

    close(FOLDER);
}

sub write
{
    my ($self) = @_;

    print "$self->{envelop}";

    foreach my $line (@{$self->{headlist}}) {
      print "$line";
    }
    print "\n";

    foreach my $line (@{$self->{body}}) {
      print "$line";
    }
}

##############################################################

sub make_hash_header
{
    my ($self,$header) = @_;

    my %hash = ();
    my @keys = ();
    my $key = "";

    my @lines = ();

    foreach my $saveline (@$header) {
      my $line = $saveline;
      push(@lines,$line);
      chomp($line);
      if( $line =~ /^(\S+):\s*(.*)$/ ) {
        $key = $1;
        if( $hash{$key} ) {	# not unique
	  #print STDERR "key is not unique: $line";
	  $hash{$key} .= "\n , " . $2;
	} else {
          $hash{$key} = $2;
	  push(@keys,$key);
	}
      } elsif( $line =~ /^\s+\S+/ ) {
        $hash{$key} .= "\n" . $line;
      } else {
        my $aux = join("",@lines);
	#print STDERR "---------- header ---------------\n";
	#print STDERR "$aux";
	#print STDERR "---------- header ---------------\n";
	print STDERR "*** error in header ... ignoring : |$line|\n";
        #die "Cannot parse header: |$line|\n";
      }
    }

    $self->{headlist} = $header;
    $self->{headhash} = \%hash;
    $self->{headkeys} = \@keys;
}

sub split_addresses
{
    my ($self,$adr) = @_;

    my @adr = ();
    my $new = "";

    $adr =~ tr/\n/ /;

    while( $adr =~ /^(.*?)([\(\<\",;])(.*)$/s ) {
      $new .= $1;
      my $match = $2;
      my $rest = $3;
      if( $match =~ /[,;]/ ) {
	if( $new = trim_line($new) ) {
	  push(@adr,$new);
	  #print STDERR "new: $new\n";
	}
	$new = "";
	$adr = $rest;
      } else {
	my ($start,$end) = get_matching($match,$rest);
	$new .= $match . $start;
	$adr = $end;
      }
    }

    $new .= $adr;
    if( $new = trim_line($new) ) {
      push(@adr,$new);
      #print STDERR "new: $new\n";
    }

    return @adr;
}

sub get_matching
{
    my ($match,$rest) = @_;

    my $search = $match;

    $search =~ tr/\{\[\(\<\"/\}\]\)\>\"/;
    $search = quotemeta($search);

    #print STDERR "matching: $match - $search\n";

    if( $rest =~ /^(.*?$search)(.*)$/s ) {
      return ($1,$2);
    } else {
      print STDERR "*** Cannot find matching: $match - $search - $rest\n";
      return ($rest,"");
    }
}

sub trim_line
{
    my ($line) = @_;

    $line =~ s/^\s+//;
    $line =~ s/\s+$//;

    return $line;
}

sub decode_mime
{
    my ($ctype) = @_;

    my ($type);
    my $rest = $ctype;

    $rest =~ s/\n/ /g;
    if( $rest =~ /^\s*(\S+\/\S+)\s*;?\s*/ ) {
      $type = $1;
      $rest = $';
    } else {
      die "Cannot parse mime type: $ctype\n";
    }

    # ggu : here decode boundary
}

##############################################################

sub is_mime
{
    my ($self) = @_;

    my $mime = $self->get_header("MIME-Version");
    $mime = trim_line($mime);
    $mime =~ s/\(.*\)//;				#strip comments

    return 0 unless( $mime =~ /^1\.0$/ );		#not a mime message

    my $ctype = $self->get_header("Content-type");
    if( not $ctype or $ctype =~ /^\s*$/ ) {
      $ctype = $self->get_header("Content-Type");
    }

    $self->{is_mime} = 1;
    $self->{mime_type} = $ctype;
    $self->{mime_boundary} = $ctype;

    decode_mime($ctype);

    return 1;
}

sub get_mime_boundary
{
    my ($self) = @_;

    return $self->{mime_boundary};
}

###################################
#&test_str(@ARGV);
###################################
1;
###################################

