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
# bibtex utilities
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use bibtex;
# 
# my $bib = new bibtex;
#
# $bib->read_file("my.bib");
#
##############################################################
#
# 25.10.2013    1.0	ggu     copied from elab.pl
# 04.11.2013    1.1	ggu     write items from external source
#
# version 1.1
#
##############################################################

use strict;

package bibtex;

##############################################################
#
# accepted types (gbib-types):
#
# JOURNAL JOURNAL_NONISI 			article
# INCOLLECTION INCOLLECTION_NONISI 		incollection
# BOOK						book
# PHDTHESIS					phdthesis
# MSTHESIS					mastersthesis
# PROC ABSTRACT					inproceedings
# TR						techreport
# MISC						misc
# 
##############################################################

sub new
{
    my $self;

    $self =     {
                         verbose        =>      0
                        ,nitems         =>      0
                        ,sitems         =>      []
                        ,items          =>      {}
                        ,bib2gbib       =>      {}
                };

    bless $self;
    $self->init_arrays();
    return $self;
}

sub init_arrays
{
    my $self = shift;

    my %bib2gbib = (
		 "article"		=>	"JOURNAL"
		,"incollection"		=>	"INCOLLECTION"
		,"inproceedings"	=>	"PROC"
		,"techreport"		=>	"TR"
		,"mastersthesis"	=>	"MSTHESIS"
		,"phdthesis"		=>	"PHDTHESIS"
		,"book"			=>	"BOOK"
	    );

    $self->{bib2gbib} = \%bib2gbib;
}

###########################################################################

sub read_file
{
    my ($self,$file) = @_;

    if ($file) {
       open(INPUT, '<', $file) or die "Cannot open $file for input\n";
    } else {
       *INPUT = *STDIN;
    }

    $self->parse_bib_items();

    close(INPUT);
}

sub write_file
{
    my ($self,$file,$items) = @_;

    if ($file) {
       open(OUTPUT, '>', $file) or die "Cannot open $file for output\n";
    } else {
       *OUTPUT = *STDOUT;
    }

    $self->write_bib_items($items);

    close(OUTPUT);
}

###########################################################################

sub surname_first {

  my $authors = shift;

  my %new = ();
  my $newname = "";

  foreach my $name (keys %$authors) {
    if( $name =~ /^(.*\.)\s+(.*)$/ ) {
      $newname = "$2 $1";
    } else {
      die "*** Cannot invert name: $name\n";
    }
    my $count = $authors->{$name};
    $new{$newname} = $count;
  }

  return \%new;
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

sub strip
{
	my $text = shift;

	$text =~ s/^\s*//;
	$text =~ s/\s*$//;

	return $text;
}

sub keep_uppercase
{
	my $value = shift;

	foreach my $upper (@::keepupper) {
	  $value =~ s/\{$upper\}/$upper/g;
	  $value =~ s/\b$upper\b/{$upper}/g;
	}

	return $value;
}

sub split_line
{
  my ($line,$nl,$gsep) = @_;	#nl is max chars per line (0 for one word/line)

  $gsep = "" unless $gsep;	#extra seperatore at beginning of line

  my @f = split(/\s+/,$line);
  my $ll = 0;
  my $new = "";
  my $sep = "";

  foreach my $item (@f) {
    my $l = length($item);
    if( $ll+$l > $nl ) {
      $sep = "\n$gsep";
      $ll = 0;
    } else {
      $sep = " ";
      $ll++;
    }
    $new .= $sep . $item;
    $ll += $l;
  }
  $new .= "\n";
  $new =~ s/^\s+//;

  return $new;
}

sub quote_entry
{
	my $text = shift;

	return "" unless $text;

	$text =~ s/([^\{])\\\"(\w)([^\}])/$1\{\\\"$2\}$3/g;

#	changes " into `` or ''

	my $in = 0;
	while( $text =~ /\"/ and not $text =~ /\\\"/ ) {
	  if( $in % 2 == 0 ) {
	    $text =~ s/\"/``/;
	  } else {
	    $text =~ s/\"/''/;
	  }
	  $in++;
	}

	if( $in % 2 == 1 ) {
	  die "*** Odd number of quotation marks found: $text\n";
	}

	return $text;
}

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

sub print_tex_authors
{
	my ($authors,$author_type) = @_;

	return "" unless $authors;

	my $line = "";

	my @f = split(/\s+and\s+/,$authors);
	my $first = shift(@f);
	my $last = pop(@f);
	$first = invert_author($first) if $author_type > 1;
	$last = invert_author($last) if $author_type == 3;
	$line = "$first";
	foreach my $name (@f) {
	  $name = invert_author($name) if $author_type == 3;
	  $line .= ", $name";
	}
	$line .= " and $last" if $last;
	$line .= " ";

	return $line;
}

sub invert_author
{
	my $name = shift;

	return "" unless $name;

	my @f = split(/\s+/,$name);
	my $lastname = pop(@f);
	if( $lastname =~ /\}$/ ) {		#composit last name
	  while( my $next = pop(@f) ) {
	    $lastname = "$next $lastname";
	    last if( $next =~ /^\{/ );
	  }
	}
	unshift(@f,"$lastname,");

	return join(" ",@f);
}

sub how_many_authors
{
	my $name = shift;

	return 0 unless $name;

	my @f = split(/\s+and\s+/,$name);
	my $n = @f;

	return $n;
}

#-------------------------------------------------------------------

sub clean_authors
{
	my $self = shift;

	my $items = $self->{sitems};
	my $n = 0;

	foreach my $item (@$items) {
	  $n++;
	  $self->clean_author($item,$n);
	}
}

sub clean_author
{
	my ($self,$item,$n) = @_;

	my @a = ();
	my @f = ();
	my $authors = $item->{author};

	@f = split(/\s+and\s+/,$authors);
	my $nf = @f;
	if( $nf > 2 ) {
	  print STDERR "$n (BIB $nf): $authors\n";
	  return;				# already in bibtex format
	} elsif( $nf == 2 ) {
	  my @ff = split(/\s*,\s*/,$f[0]);	# see if not final and
	  my $nff = @ff;
	  if( $nff <= 1 ) {			# just one name - ok
	    print STDERR "$n (BIB 2): $authors\n";
	    return;				# already in bibtex format
	  } else {				# is final and -> convert
	    print STDERR "$n (ORIG $nf $nff): $authors\n";
	    $authors =~ s/\s*and/,/;
	    print STDERR "$n (CONVERT $nf $nff): $authors\n";
	  }
	}

	@f = split(/\s*,\s*/,$authors);
	print STDERR "$n: $authors\n";

	while( my $item = shift(@f) ) {
	}
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

sub write_bib_items
{
	my ($self,$items) = @_;

	$items = $self->{sitems} unless $items;
	my $n = 0;

	foreach my $item (@$items) {
	  $n++;
	  my $line = $self->write_bib_item($item,$n);
	  print OUTPUT "$line\n";
	}

	print STDERR "$n items written\n";
}

sub write_bib_item
{
	my ($self,$item,$n) = @_;

	my $keys = $self->check_item($item,$n);

	my $name = $item->{"ggu_bibtype"};
	$::number++;
	my $line = "\@$name\{bibentry_$::number\n";

	foreach my $key (@$keys) {
	  $line .= write_bib_key($item,$key);
	}

	$line .= "\t}\n";

	return $line;
}

sub write_bib_key
{
	my ($entry,$key) = @_;

	my $value = $entry->{$key};
	if( $::keep_uppercase and $key eq "title" ) {
	  #$value = keep_uppercase($value);
	}
	if( $value ) {
	  $value = split_line($value,50,"\t\t\t ");
	  $value =~ s/\n$//;
	  my $line = "\t,$key = \t\"$value\"\n";
	  return $line;
	} else {
	  return "";
	}
}

sub check_item
{
        my ($self,$item,$n) = @_;

        my @keys = qw/
                        author year title pages 
                        journal volume number
                        booktitle editor publisher series
                        type institution
                        school
                        note rest
			issn isbn doi url
                        ggu_type ggu_link
                        /;

        my %entry_count = ();

        $entry_count{"ggu_content"}++;
        $entry_count{"ggu_bibtype"}++;
        $entry_count{"ggu_bibkey"}++;

        $entry_count{"affiliation"}++;
        $entry_count{"keywords"}++;
        $entry_count{"address"}++;
        $entry_count{"language"}++;
        $entry_count{"abstract"}++;

        #$entry_count{"ggu_type"}++;
        #$entry_count{"ggu_link"}++;

        foreach my $key (@keys) {
          $entry_count{$key}++;
        }

        foreach my $key (keys %$item) {
          unless( $entry_count{$key} ) {
	    if( $key =~ /^CNR_/ ) {
	      ;		#ok
	    } else {
	      my $name = $item->{ggu_bibkey};
              die "*** Key $key in entry $name ($n) not used...\n";
	    }
          }
        }

        return \@keys;
}

#-------------------------------------------------------------------

sub parse_bib_items
{
	my $self = shift;

	my $in_item = 0;
	my $line = "";

	#print "......... parsing bib file...\n";

	while(<INPUT>) {
	  chomp;
	  if( /^\s*@/ ) {		# new item
	    if( $in_item ) {
		$self->parse_bib_item($line);
		$in_item = 0;
		$line = "";
	    }
	    #print "new item found: $_\n";
	    $line = "$_ ";
	    $in_item = 1;
	  } else {	
	    $line .= " $_"
	  }
	}

	$self->parse_bib_item($line) if $in_item;

	my $nitems = $self->{nitems};
	print STDERR "$nitems items read\n";
}

sub parse_bib_items_0
{
	my $self = shift;

	my $in_item = 0;
	my $line = "";

	#print "......... parsing bib file...\n";

	while(<INPUT>) {
	  chomp;
	  if( $in_item ) {
	    if( /^\s*$/ ) {		# empty line -> finshed scanning item
	        #print "new item finished\n";
		$self->parse_bib_item($line);
		$in_item = 0;
		$line = "";
	    } else {
	      $line .= " $_"
	    }
	  } elsif( /^\s*@/ ) {		# new item
	    if( $in_item ) {
	      die "already in item: cannot open another one: $_\n";
	    }
	    #print "new item found: $_\n";
	    $line = "$_ ";
	    $in_item = 1;
	  } else {			# comment
	    ;
	  }
	}

	$self->parse_bib_item($line) if $in_item;

	my $nitems = $self->{nitems};
	print STDERR "$nitems items read\n";
}

sub parse_bib_item
{
	my ($self,$line) = @_;

	my $type;
	my $bibkey;
	my $rest;
	my $key;
	my $value;

	my $orig = $line;

	#print "new bib item:\n";
	#print "$line\n";

	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	$line =~ s/\s+/ /g;

	if( $line =~ s/\s*\}$// ) {		#chop ending }
	  ;
	} else {
	  die "Cannot parse bib item ending: $line\n";
	}

	if( $line =~ /^@(\w+)\s*\{\s*(\w+)\s*/ ) {	#get header
	  $type = $1;
	  $bibkey = $2;
	  $rest = $';
	} else {
	  die "Cannot parse bib item header: $line\n";
	}

	my $item = $self->init_entry($type,$bibkey);

	#print STDERR "bib item: $type   $bibkey\n";

	while( $rest ) {
	  if( $rest =~ /^\s*,\s*(\w+)\s*=\s*/ ) {
	    $key = $1;
	    $rest = $';
	    ($value,$rest) = get_pars($rest);
	    #print "key: $key = $value\n";
	    $value = pop_pars($value);
	    $self->add_entry($item,$key,$value);
	  } elsif( $rest =~ /\s*,\s*/ ) {		# trailing comma
	    $rest = "";
	  } else {
	    print STDERR "$orig\n";
	    die "Cannot parse bib item keys: $rest\n";
	  }
	}

	$self->insert_entry($item);
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub init_entry
{
	my ($self,$type,$bibkey) = @_;

	my %item = ();
	$item{"ggu_bibtype"} = $type;
	$item{"ggu_bibkey"} = $bibkey;

	return\%item;
}

sub add_entry
{
        my ($self,$item,$key,$value) = @_;

        $value = quote_entry($value);

        $item->{$key} = $value if $value;
}

sub insert_entry
{
        my ($self,$item) = @_;

	my $sitems = $self->{sitems};
	push(@$sitems,$item);

	my $items = $self->{items};
	my $bibkey = $item->{"ggu_bibkey"};
	$items->{$bibkey} = $item;

	$self->{nitems}++;
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub get_pars
{
	my $line = shift;

	my $open = substr($line,0,1);
	my $close;

	if( $open eq '"' ) {
	  $close = $open;
	} elsif( $open eq '{' ) {
	  $close = '}';
	} else {
	  die "cannot handle parenthesis: $open\n";
	}

	# here we still have to handle things like \" \{ \}

	my $content = $open;
	my $nc = 1;
	my $rest = substr($line,1);	# without opening
	while( $rest and $nc ) {
	  if( $rest =~ /^(.*?)([\"\{\}])/ ) {
	    my $text = $1;
	    my $before = substr($text,-1,1);
	    my $what = $2;
	    if( $before eq '\\' ) {
	      ;
	    } elsif( $what eq $close ) {
	      $nc--;
	    } elsif( $what eq $open ) {
	      $nc++;
	    }
	    $content .= $1 . $what;
	    $rest = $';
	  } else {			# read until end of line
	    $content .= $rest;
	    $rest = "";
	  }
	}

	if( $nc ) {
	  die "cannot find closing pars: $line\n";
	} 

	return ($content,$rest);
}

sub pop_pars
{
	my $value = shift;

	$value =~ s/^[\{\}\"]//;
	$value =~ s/[\{\}\"]$//;

	return $value;
}

#--------------------------------------------------------------
1;
#--------------------------------------------------------------

