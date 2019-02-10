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
# utilities for reading and writing html files
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use html3;
# 
# my $html = new html3;
# $html->open_file_write("g.html");	
# $html->write_header("title");
# $html->write_trailer();
# $html->close_file_write();
#
# my $html = new html3;
# $html->read_file("g.html");	
# my $rtext = $html->{text};
# my $text = $$rtext;
# my ($content,$options,$rest) = $html->get_tag("tr",$text);
#
#---------------------------------------
#
# version 2.0	14.10.2010	new routines
# version 3.0	12.03.2016	new read routines
# version 3.1	06.12.2018	documentation and some more routines
#
##############################################################

use strict;

package html3;

##############################################################

sub new
{
    my ($pck) = @_;

    my $self;

    $self =	{
	    		  file_name_read	=>	undef
	    		 ,file_handle_read	=>	undef
	    		 ,file_name_write	=>	undef
	    		 ,file_handle_write	=>	undef
			 ,text			=>	undef
			 ,list			=>	[]
			 ,version		=>	"3.1"
		};

    bless $self;
    return $self;
}

#-------------------------------------------------------------------
# file opening routines
#-------------------------------------------------------------------

sub open_file_write
{
    my ($self,$file) = @_;

    if( $file ) {
      open(FILE,">$file") or die "Cannot open file: $file\n";
      $self->{file_handle_write} = \*FILE;
      $self->{file_name_write} = $file;
    } else {
      $self->{file_handle_write} = \*STDOUT;
      $self->{file_name_write} = "-";
    }
}

sub open_file_read
{
    my ($self,$file) = @_;

    if( $file ) {
      open(FILE,"<$file") or die "Cannot open file: $file\n";
      $self->{file_handle_read} = \*FILE;
      $self->{file_name_read} = $file;
    } else {
      $self->{file_handle_read} = \*STDIN;
      $self->{file_name_read} = "-";
    }
}

sub read_file
{
    my ($self,$file) = @_;

    $self->open_file_read($file);
    local $/=undef;
    my $FH = $self->{file_handle_read};
    my $text = <$FH>;
    $self->{text} = \$text;
    $self->close_file_read($file);
}

sub close_file_write
{
    my ($self) = @_;

    my $file_handle = $self->{file_handle_write};
    my $file_name   = $self->{file_name_write};

    if( $file_name ne "-" ) {
      close($file_handle);
      $self->{file_handle} = undef;
      $self->{file_name} = undef;
    }
}

sub close_file_read
{
    my ($self) = @_;

    my $file_handle = $self->{file_handle_read};
    my $file_name   = $self->{file_name_read};

    if( $file_name ne "-" ) {
      close($file_handle);
      $self->{file_handle} = undef;
      $self->{file_name} = undef;
    }
}

sub print_version {

  my ($self) = @_;
  my $version = $self->{version};
  print STDERR "html.pm version $version\n";
}

#-------------------------------------------------------------------
# writing file
#-------------------------------------------------------------------

sub write_header {

  my ($self,$title) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<html>\n";
  print $FH "<head>\n";
  print $FH "<title>$title</title>\n";
  print $FH "</head>\n";
  print $FH "<body>\n";
  print $FH "\n";
}

sub write_trailer {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "\n";
  print $FH "</body>\n";
  print $FH "</html>\n";
}

sub write_header_redirect {

  my ($self,$title,$html,$seconds) = @_;

# use $html as: url=http://www.yourdomain.com/index.html

  $seconds = 0 unless $seconds;

  my $FH = $self->{file_handle_write};

  print $FH "<html>\n";
  print $FH "<head>\n";
  print $FH "<title>$title</title>\n";
  print $FH "<meta HTTP-EQUIV=\"REFRESH\" content=\"$seconds; url=$html\">";
  print $FH "</head>\n";
  print $FH "<body>\n";
  print $FH "\n";
}

#------------------------------------------------------------------- table

sub open_table {

  my ($self,$extra) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "\n";
  if( $extra ) {
    print $FH "<table $extra>\n";
  } else {
    print $FH "<table>\n";
  }
}

sub close_table {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "</table>\n";
  print $FH "\n";
}

sub open_table_row {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<tr>\n";
}

sub close_table_row {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "</tr>\n";
}

sub insert_table_data {

  my ($self,$data,$color) = @_;

  my $FH = $self->{file_handle_write};

  if( $color ) {
    print $FH "<td><font color=\"$color\">$data</font></td>\n";
  } else {
    print $FH "<td>$data</td>\n";
  }
}

#-------------------------------------------------------------------

sub open_list {

  my ($self,$type,$text) = @_;

  # type can be dir/menu/ol/ul - ul is default - can have options

  my $FH = $self->{file_handle_write};

  $type = "ul" unless $type;

  print $FH "\n";
  print $FH "<$type>\n";
  print $FH "<hl>$text</hl>\n" if $text;	# header

  $type =~ s/\s+.*$//;
  my $list = $self->{list};
  push(@$list,$type);
}

sub close_list {

  my ($self) = @_;

  my $list = $self->{list};
  my $type = pop(@$list);

  my $FH = $self->{file_handle_write};

  print $FH "</$type>\n";
  print $FH "\n";
}

sub insert_list {

  my ($self,$text) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<li>$text</li>\n";
}

#-------------------------------------------------------------------

sub insert_heading {

  my ($self,$data,$level) = @_;

  my $FH = $self->{file_handle_write};

  $level = 1 unless $level;
  my $tag = "h" . $level;
  print $FH "<$tag>$data</$tag>\n";
}

sub insert_para {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<p>$data</p>\n";
}

sub insert_data {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "$data";
}

sub insert_line {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "$data\n";
}

sub insert_break {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<br>\n";
}

sub insert_image {

  my ($self,$image,$options) = @_;

  my $FH = $self->{file_handle_write};

  $options = "" unless $options;

  print $FH "<img src=\"$image\" $options>";
}

#-------------------------------------------------------------------

sub make_anchor {

  my ($self,$text,$href) = @_;

  my $line = "<a href=\"$href\">$text</a>";

  return $line;
}

sub make_clickable_image {

  my ($self,$img_file,$options,$href) = @_;

  my $img_line = "<img src=\"$img_file\" $options>";
  my $line = "<a href=\"$href\">$img_line</a>";

  return $line;
}

#-------------------------------------------------------------------
# reading html tags
#-------------------------------------------------------------------

sub find_next_tag { # finds and returns next available tag without reading it

  my ($self,$text) = @_;

  if( $text =~ /^(.*?)<(\w+)(\W)/s ) {
    my $tag = $2;
    return $tag;
  } else {
    return;
  }
}

sub get_nth_tag {

# gets nth specified tag, else as get_tag (discards n-1 tags)

  my ($self,$tag,$text,$n) = @_;

  $n = 0 unless $n;
  my ($content,$options);

  do {
    ($content,$options,$text) = $self->get_tag($tag,$text);
    $n = 0 unless $content;
  } while( --$n > 0 );

  return ($content,$options,$text);
}

sub get_tag { 

# gets next specified tag in text and returns content, options of tag and rest

  my ($self,$tag,$text) = @_;

  my $closing = "";
  my $before = "";
  my $options = "";
  my $contents = "";

  #my $line = substr($text,0,30);
  #print STDERR "reading this: $line\n";

  # find opening tag ----------------------------------

  if( $text =~ /^(.*?)<$tag(\W)/si ) {
    $before = substr($1,-30);
    $closing = $2;
    $text = $';
  } else {
    return
  }
  #my $after = substr($text,0,20);
  #print STDERR "have found... |$closing|$before|$after|\n";

  # find options in tag ----------------------------------

  if( $closing eq ">" ) {		# full tag read, no options
    # wait to read contents
    #print STDERR "have found without options...\n";
  } elsif( $closing eq " " ) {	# space read - must still read possible options
    #print STDERR "have found with options...\n";
    if( $text =~ /^(.*?)>/s ) {
      $options = $1;
      $text = $';
    } else {
      die "cannot find closing > for tag $tag\n";
    }
  } else {
    my $line = substr($text,0,20);
    die "Error reading tag $tag |$closing|: $line...\n";
  }
    
  # find closing tag ----------------------------------

  if( $text =~ /^(.*?)(<\/$tag>)/si ) {
    $contents = $1;
    my $matched = $2;
    my $line1 = substr($1.$matched,-30);
    #print STDERR "have found end tag... $line1\n";
    $text = $';
    my $line2 = substr($text,0,30);
    #print STDERR "next text... $line2\n";
    #print STDERR "--------------------------\n";
  } else {
    my $line = substr($text,0,20);
    die "Error reading end tag $tag: $line...\n";
  }

  return ($contents,$options,$text);
}

sub get_content { # gets content of tag, throws away rest

  my ($self,$tag,$text) = @_;

  my ($content,$options,$rest) = $self->get_tag($tag,$text);

  return $content;
}

sub delete_tag { # not yet working

  my ($self,$text) = @_;
  my $tag;

  my ($content,$options,$rest) = $self->get_tag($tag,$text);

  return $content;
}

sub clean_tag {	# compress white space

  my ($self,$text) = @_;

  return "" unless $text;

  $text =~ s/\n/ /sg;
  $text =~ s/^\s+//;
  $text =~ s/\s+$//;
  $text =~ s/\s+/ /sg;

  return $text;
}

sub show_text { # shows limited number of chars of text

  my ($self,$text,$n) = @_;

  $n = 40 unless $n;
  my $ptext = substr($text,0,$n);

  print "$ptext\n";
}

#-------------------------------------------------------------------

################################
1;
################################
