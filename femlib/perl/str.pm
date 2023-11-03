#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# str utility routines
#
# version	1.10    02.11.2023      die if section is not found
# version	1.9     30.10.2023      bug fix for writing section with decrp
# version	1.8     26.06.2023      new routine has_section()
# version	1.7     11.06.2023      parameter nocomment implemented
# version	1.6     08.06.2023      write extra section with description
# version	1.5     30.03.2023      prepared new number section
# version 	1.4	16.01.2022	parses extra section with description
# version 	1.3	24.11.2014	finds section now, can insert values
# version 	1.2	11.05.2011	restructured and commented
# version 	1.1	?
#
# Usage:
#
# #!/usr/bin/env perl
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
# use str;
# 
# my $file = $ARGV[0];
# 
# my $str = new str;
# $str->read_str($file);
# $str->print_sections();
# $str->write_str();
#
#------------------------------------------------
#
# todo:
#
# write number section better (more compact)
# finish write_line_80() and use while writing array, number, param section
#
##############################################################

use strict;

package str;

##############################################################

%str::namesections = ( 
		 para	=>	1
		,bound	=>	1
		);

##############################################################

sub new
{
    my $self;

    $self =	{
	    		 file		=>	undef
			,sections	=>	{}
			,sequence	=>	[]
			,section_types	=>	undef
			,verbose	=>	0
			,quiet		=>	0
			,nocomment	=>	0
		};

    bless $self;
    $self->make_arrays();
    return $self;
}

#-----------------------------------------------------------------
# internal routines
#-----------------------------------------------------------------


sub make_arrays {

  my ($self) = @_;

  my @param_sections = qw/ para bound name color legvar wrt sedtr waves close /;
  my @number_sections = qw/ extra flux levels /;
  my @table_sections = qw/ area extts /;
  my @title_sections = qw/ title /;
  my @comment_sections = qw/ comment /;

  my %sections = ();

  $self->insert_arrays(\%sections,"param",\@param_sections);
  $self->insert_arrays(\%sections,"number",\@number_sections);
  $self->insert_arrays(\%sections,"table",\@table_sections);
  $self->insert_arrays(\%sections,"title",\@title_sections);
  $self->insert_arrays(\%sections,"comment",\@comment_sections);

  $self->{section_types} = \%sections;
}

sub insert_arrays {

  my ($self,$sections,$type,$array) = @_;

  foreach my $section (@$array) {
    $sections->{$section} = $type;
  }
}

sub delete_str {

  my ($self) = @_;

  $self->{file} = undef;
  $self->{sections} = {};
  $self->{sequence} = [];
}

#-----------------------------------------------------------------
# accessor routines
#-----------------------------------------------------------------

sub set_verbose {

  my ($self,$value) = @_;

  $self->{verbose} = $value;
}

sub set_nocomment {

  my ($self,$value) = @_;

  $self->{nocomment} = $value;
}

# section defaults to "para"
# number defaults to undef

sub get_title {  return $_[0]->{title}; }
sub get_simul {  return $_[0]->{simul}; }
sub get_basin {  return $_[0]->{basin}; }

sub set_title {  $_[0]->{title} = $_[1]; }
sub set_simul {  $_[0]->{simul} = $_[1]; }
sub set_basin {  $_[0]->{basin} = $_[1]; }

sub get_numbers {

  my ($self,$section,$number) = @_;

  my $sect = $self->get_section($section,$number);

  return $sect->{array};
}

sub set_numbers {

  my ($self,$numbers,$section,$number) = @_;

  my $sect = $self->get_section($section,$number);

  $sect->{array} = $numbers;
}

sub add_numbers {

  my ($self,$numbers,$section,$number) = @_;

  my $sect = $self->get_section($section,$number);

  my $array = $sect->{array};
  if( ref($numbers) eq "ARRAY" ) {
    push(@$array,@$numbers);
  } else {
    push(@$array,$numbers);
  }
}

sub get_value {

  my ($self,$name,$section,$number) = @_;

  my $sect = $self->get_section($section,$number);
  my $hash = $sect->{hash};

  return $hash->{$name};
}

sub set_value {

  my ($self,$name,$value,$section,$number) = @_;

  my $sect = $self->get_section($section,$number);
  my $hash = $sect->{hash};

  if( not defined $hash->{$name} ) {
    my $items = $sect->{items};
    push(@$items,$name);
  }
  $hash->{$name} = $value;
}

sub get_section {

  my ($self,$section_name,$number) = @_;

  my $section_id = "para";	#section defaults to this
  $section_id = $section_name if defined $section_name;
  $section_id .= $number if defined $number;

  my $sections = $self->{sections};
  my $sequence = $self->{sequence};

  foreach my $section (@$sequence) {
	  my $sect = $sections->{$section};
	  return $sect if( $sect->{id} eq $section_id );
  }

  die "*** Cannot find section: $section_id\n";

  return;
}

#-----------------------------------------------------------------
# prepare sections and comments
#-----------------------------------------------------------------

sub make_section {

  my ($self,$name,$info) = @_;

  my $number = "";
  if( $name =~ /([a-zA-Z_]+)(\d+)/ ) {
	  $name = $1;
	  $number = $2;
  }
  my $section_id = "$name$number";

  my $sect =    {
                         id             =>      $section_id
                        ,name           =>      $name
                        ,number         =>      $number
                        ,info           =>      $info
                        ,data           =>      []	# original lines
                        ,hash           =>      {}	# hash of items
                        ,items          =>      []	# names of items
                };

  my $sections = $self->{sections};
  $sections->{$section_id} = $sect;

  my $sequence = $self->{sequence};
  push(@$sequence,$section_id);

  my $section_types = $self->{section_types};
  my $type = $section_types->{$name};

  print STDERR "section: $name ($type) - $number - $info\n" if $self->{verbose};

  return $sect;
}

sub make_comment {

  my ($self,$number,@data) = @_;

  $number++;
  my $name = "comment$number";

  my $sect = $self->make_section($name,"");
  $sect->{data} = \@data;

  return $number;
}

#-----------------------------------------------------------------
# read str and sections
#-----------------------------------------------------------------

sub read_str {

  my ($self,$file) = @_;

  my $comment = 0;
  my @comment = ();

  return unless $file;
  open(STR_FILE,"$file") || die "Cannot open file: $file\n";

  $self->delete_str();	#erase info on last str file

  $self->{file} = $file;
  print STDERR "reading file: $file\n" if $self->{verbose};

  while( <STR_FILE> ) {

	  chomp;
	  s/\r$//;

	  if( /^\s*[\$\&](\w+)\s*(.*)/ ) {

		  $comment = $self->make_comment($comment,@comment);
		  @comment = ();

		  my $section = $1;
		  my $section_info = $2;

		  my $sect = $self->make_section($section,$section_info);

		  $sect->{data} = $self->read_section($sect);

		  $self->parse_section($sect);		#must still insert

	  } else {	#comment
		chomp;
		push(@comment,$_);
	  }

  }

  $comment = $self->make_comment($comment,@comment);

  close(STR_FILE);
}

sub read_section {

  my ($self,$sect) = @_;

  my $section = $sect->{name};
  my $number = $sect->{number};
  $number = "" unless $number;
  $number = "($number)" if $number;
  print STDERR "reading section $section $number\n" if $self->{verbose};;

  my @data = ();

  while(<STR_FILE>) {
	  chomp;
	  if( /^\s*[\$\&]end\s*$/ ) {
		  return \@data;
	  } elsif( /^\s*[\$\&](\w+)\s*(.*)/ ) {
		  die "new section $1 started before ending $section\n";
	  }
	  push(@data,$_);
  }
  die "end of section $section not found\n";
}

#-----------------------------------------------------------------
# parse sections
#-----------------------------------------------------------------

sub parse_section {

  my ($self,$sect) = @_;

  my $section_name = $sect->{name};
  my $section_types = $self->{section_types};
  my $type = $section_types->{$section_name};

  if( $type ) {
    #print STDERR "The section is of type $type\n";
  } else {
    print STDERR "*** Unknown type of section: $section_name\n";
    $type = "unknown";
    #die "*** Unknown type of section: $section_name\n";
  }

  if( $type eq "param" ) {
    $self->parse_param_section($sect);
  } elsif( $type eq "title" ) {
    $self->parse_title_section($sect);
  } elsif( $type eq "number" ) {
    $self->parse_number_section($sect);
  } elsif( $type eq "table" ) {
    $self->parse_table_section($sect);
  } else {
    $self->skip_section($sect);
    #die "*** Cannot yet parse section type $type\n";
  }
}

sub parse_param_section {

  my ($self,$sect) = @_;

  my %hash = ();
  my $item;
  my $debug = 0;
  my $data = $sect->{data};
  my $items = $sect->{items};
  my $name = $sect->{name};
  my $dline = join(" ",@$data);
  $dline .= "  ";	# just to be sure that there is some ws at the end

  while( $dline =~ /^\s*(\w+)\s*=\s*/ ) {
	my $name = $1;
	my $rest = $';
	my @val = ();
	while( defined ($item = $self->get_next_value(\$rest)) ) {
	  push(@val,$item);
	}
	my $n = @val;
	if( $n == 1 ) {		#only one value
	  $hash{$name} = $val[0];
	  push(@$items,$name);
	  print "new value: $name  $val[0]\n" if $debug;
	} elsif( $n > 1 ) {	#array
	  my $line = join(" ",@val);
	  print "new array: $name  $line\n" if $debug;
	  $hash{$name} = \@val;
	  push(@$items,$name);
	} else {
	  die "Error in parsing STR file: for $name $n values found\n";
	}
	$dline = $rest;
  }

  $sect->{hash} = \%hash;
}

sub get_next_value {

  my ($self,$rline) = @_;

  my $val = undef;
  my $first = substr($$rline,0,1);

  if( $$rline =~ /^\s*$/ ) {			# nothing more
  } elsif( $first =~ /[a-z_]/i ) {		# new name
  } elsif( $first eq "\'" ) {			# string
    if( $$rline =~ /^\'(.*?)\'[\s,]\s*/ ) {
      $val = "'$1'";
      $$rline = $';
    }
  } elsif( $$rline =~ /^([0-9+-\.ed]+)\s*/i ) {	# number
    $val = $1;
    $$rline = $';
  } else {
    die "Cannot parse: $$rline\n";
  }

  if( $$rline =~ /^,\s*/ ) {			# strip comma
    $$rline = $';
  }

  return $val;
}

sub parse_title_section {

  my ($self,$sect) = @_;

  my $data = $sect->{data};

  $self->{title} = $data->[0];
  $data->[1] =~ s/^\s+//;
  $self->{simul} = $data->[1];
  $data->[2] =~ s/^\s+//;
  $self->{basin} = $data->[2];
}

sub parse_number_section {

  my ($self,$sect) = @_;

  my ($dline,$description);

  my $data = $sect->{data};
  my $name = $sect->{name};
  my @number = ();
  my @description = ();

  print STDERR "parsing number section... $name\n" if $self->{verbose};

  foreach my $line (@$data) {
    if( $line =~ /(.+)\'(.+)\'/ ) {	# with description... ' ' must be last
      $dline = $1;
      $description = $2;
    } else {
      $dline = $line;			# do not touch data section
      $description = "";
    }
    $dline =~ s/,/ /g;		# subst comma with white space
    $dline =~ s/^\s+//;
    $dline =~ s/\s+$//;
    my @f = split(/\s+/,$dline);
    push(@number,@f);
    push(@description,$description) if $description;
    if( $description and $name eq "flux" ) {
      push(@number,0);		#insert 0 to divide section
    }
  }

  my $n = @number;
  my $d = @description;
  my $b = count_blocks(\@number);
  if( $d > 0 ) {
    if( $name eq "extra" and $n != $d ) {
      print STDERR "number of descriptions not matching numbers given\n";
      die "error in section $name: $d $n\n";
    } elsif( $name eq "flux" and $b != $d ) {
      print STDERR "number of descriptions not matching blocks given\n";
      die "error in section $name: $d $b\n";
    } elsif( $name eq "levels" ) {
      print STDERR "no description allowed in this section\n";
      die "error in section $name: $d $n\n";
    }
  }

  if( $name eq "extra" ) {
    $sect->{array} = \@number;
    $sect->{description} = \@description;
  } elsif( $name eq "flux" ) {
    $sect->{array} = \@number;
    $sect->{description} = \@description;
  } elsif( $name eq "levels" ) {
    $sect->{array} = \@number;
    $sect->{description} = \@description;
  }

  return;

  # this is a debug section

  print STDERR "section $name: $n $b $d\n";

  foreach my $dd (@description) {
    print STDERR "$dd\n";
  }
  foreach my $nn (@number) {
    print STDERR "$nn\n";
  }
  
}

sub parse_table_section {

  my ($self,$sect) = @_;

  # no parsing
}

sub skip_section {

  my ($self,$sect) = @_;

  # no parsing
}

sub count_blocks {

  my $list = shift;

  my $nblocks = 0;
  my $in_block = 0;

  foreach my $n (@$list) {
    $in_block = 0 if $n == 0;
    if( $n > 0 and $in_block == 0 ) {
      $in_block = 1;
      $nblocks++;
    }
  }

  return $nblocks;
}

#-----------------------------------------------------------------
# apply function to sections
#-----------------------------------------------------------------

sub apply_sections {

  my ($self,$func) = @_;

  my $sections = $self->{sections};
  my $sequence = $self->{sequence};

  foreach my $section (@$sequence) {
	  my $sect = $sections->{$section};
	  &$func($self,$sect);
  }
}

#-----------------------------------------------------------------
# print sections (for debug)
#-----------------------------------------------------------------

sub print_sections {

  my ($self) = @_;

  my $sections = $self->{sections};
  my $sequence = $self->{sequence};

  foreach my $section (@$sequence) {
	  my $sect = $sections->{$section};
	  $self->print_section($sect);
  }

}

sub print_section {

  my ($self,$sect) = @_;

  my $id		=	$sect->{id};
  my $name		=	$sect->{name};
  my $number		=	$sect->{number};
  my $info		=	$sect->{info};
  my $data		=	$sect->{data};

  print "$id: $name $number  $info\n";

  foreach my $line (@$data) {
	  print "$line\n";
  }
}

#-----------------------------------------------------------------
# write str and sections
#-----------------------------------------------------------------

sub write_str {

  my ($self,$file) = @_;

  if( $file ) {
    open(STR_FILE,">$file") || die "Cannot open file $file for writing\n";
    #$oldfh = select();
    select(STR_FILE);
  }

  $self->write_sections();

  if( $file ) {
    select(STDOUT);
  }
}

sub write_sections {

  my ($self) = @_;

  my $sections = $self->{sections};
  my $sequence = $self->{sequence};

  foreach my $section (@$sequence) {
	  my $sect = $sections->{$section};
	  $self->write_section($sect);
  }

}

sub write_section {

  my ($self,$sect) = @_;

  my $id		=	$sect->{id};
  my $name		=	$sect->{name};
  my $number		=	$sect->{number};
  my $info		=	$sect->{info};
  my $data		=	$sect->{data};
  my $type		=	$self->{section_types}->{$name};

  $type = "unknown" unless $type;

  $self->write_start_of_section($sect);

  #print "+++++++++ type: $type  $name  ($data)\n";

  if( $type eq "param" ) {
    $self->write_param_section($sect);
  } elsif( $type eq "comment" ) {
    if( $self->{nocomment} == 0 ) {
      $self->write_comment_section($sect);
    }
  } elsif( $type eq "title" ) {
    $self->write_title_section($sect);
  } elsif( $type eq "number" ) {
    $self->write_number_section($sect);
  } elsif( $type eq "table" ) {
    $self->write_table_section($sect);
  } else {
    $self->write_unknown_section($sect);
    #die "*** Cannot yet write section type $type\n";
  }

  $self->write_end_of_section($sect);
}

sub write_start_of_section {

  my ($self,$sect) = @_;

  my $name		=	$sect->{name};
  my $number		=	$sect->{number};
  my $info		=	$sect->{info};
  my $type		=	$self->{section_types}->{$name};

  $type = "unknown" unless $type;

  return if $type eq "comment";
  
  my $line = "$name$number";
  $line .= "  $info" if $info;
  print "\$$line\n";
}

sub write_end_of_section {

  my ($self,$sect) = @_;

  my $name		=	$sect->{name};
  my $type		=	$self->{section_types}->{$name};

  $type = "unknown" unless $type;

  return if $type eq "comment";

  print "\$end\n";
}

sub write_comment_section {

  my ($self,$sect) = @_;

  my $data		=	$sect->{data};
  my $n = @$data;

  foreach my $line (@$data) {
	  print "$line\n";
  }
}

sub write_param_section {

  my ($self,$sect) = @_;

  my $hash = $sect->{hash};

  #my @items = sort keys %$hash;
  my @items = @{$sect->{items}};

  foreach my $name (@items) {
    my $value = $hash->{$name};
    if( ref($value) eq "ARRAY" ) {
      write_array($value,"\t$name = ");
    } else {
      print "\t$name = $value\n";
    }
  }
}

sub write_title_section {

  my ($self,$sect) = @_;

  my $hash = $sect->{hash};

  print "$self->{title}\n";
  print "\t$self->{simul}\n";
  print "\t$self->{basin}\n";
}

sub write_number_section {

  my ($self,$sect) = @_;

  my $name = $sect->{name};
  my $data = $sect->{array};
  my $description = $sect->{description};
  my $na = @$data;
  my $nd = @$description;
  my $nb = count_blocks($data);

  #print STDERR " blocks: $name $na $nd $nb\n";

  if( $nd == 0 ) {
    write_array($data);
  } elsif( $na == $nd ) {	# section without blocks
    write_array_with_description($data,$description);
  } elsif( $nb == $nd ) {	# section with blocks
    write_array_blocks_with_description($data,$description);
  } else {
    die "cannot write number section $name: $na $nd $nb\n";
  }
}

sub write_array_blocks_with_description {

  my ($array,$description) = @_;

  my $nval = 5;			# how many values on one line

  my $na = @$array;
  my $nd = @$description;

  #print STDERR "array $na\n";
  #foreach my $f (@$array) {
  #  print STDERR "$f ";
  #}
  #print STDERR "\n";
  #print STDERR "description $nd\n";
  #foreach my $f (@$description) {
  #  print STDERR "$f ";
  #}
  #print STDERR "\n";

  my $j = 0;
  my $in_block = 0;
  my $to_print = 0;
  my $line = "";

  for( my $i=0; $i<$na; $i++ ) {
    my $item = $array->[$i];
    if( $item == 0 ) {
      next if( $in_block == 0 );
      $j = 0;
      $in_block = 0;
      my $descr = shift(@$description);
      $line .= "    '" . $descr . "'";
      $to_print = 1;
    } else {
      $j++;
      $in_block = 1;
      $line .= "  $item";
      if( $j == $nval ) {
        $to_print = 1;
        $j = 0;
      }
    }
    if( $to_print ) {
      print "$line\n";
      $to_print = 0;
      $line = "";
    }
  }
}

sub write_array_with_description {

  my ($array,$description) = @_;

  my $na = @$array;

  for( my $i=0; $i<$na; $i++ ) {
    my $item = $array->[$i];
    my $descr = $description->[$i];
    $descr = "'" . $descr . "'";
    print "\t$item\t$descr\n";
  }
}

sub write_array {

  my ($array,$extra) = @_;

  my $nval = 5;			# how many values on one line
  my $line = "";

  print "$extra" if $extra;

  my $na = @$array;

  my $i = 0;
  foreach my $item (@$array) {
    $i++;
    print "   $item";
    if( $item == 0 or $i == $nval ) {
      print "\n";
      $i = 0;
    }
  }
  print "\n" unless $i == 0;
}

sub write_table_section {

  my ($self,$sect) = @_;

  my $data = $sect->{data};

  foreach my $line (@$data) {
    print "$line\n";
  }
}

sub write_unknown_section {	# simply copy

  my ($self,$sect) = @_;

  my $data = $sect->{data};

  foreach my $line (@$data) {
    print "$line\n";
  }
}

#-----------------------------------------------------------------
# utilities
#-----------------------------------------------------------------

sub has_section {

  my ($self,$sectname) = @_;

  my $sections = $self->{sections};
  my $sequence = $self->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};
    my $name = $sect->{name};
    #print STDERR "section $name\n";
    return 1 if $name eq $sectname;
  }
  return 0;
}

sub write_line_80
{
  my $line = shift;

  # not yet ready...
}

#-----------------------------------------------------------------
# testing
#-----------------------------------------------------------------

sub test_str
{
    my @files = @_;

    my $str = new str;

    print "reading...\n";
    $str->read_str($files[0]);
    print "writing...\n";
    $str->print_sections;
}

#&test_str(@ARGV);

#-----------------------------------------------------------------
# return 1;
#-----------------------------------------------------------------

1;

#-----------------------------------------------------------------
# end of routine
#-----------------------------------------------------------------

