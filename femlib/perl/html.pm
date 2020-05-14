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
# utilities for writing html files
#
# example of usage:
#
# use html;
# 
# my $html = new html;
# $html->open_file("g.html");	# also: my $html = new html("g.html);
# $html->write_header("title");
# $html->write_trailer();
#
#---------------------------------------
#
# version 2.0	14.10.2010	new routines
#
##############################################################

use strict;

package html;

##############################################################

sub new
{
    my ($pck,$file) = @_;

    my $self;

    $self =	{
	    		  file_name	=>	undef
	    		 ,file_handle	=>	undef
			 ,list		=>	[]
			 ,version	=>	"2.0"
		};

    bless $self;
    $self->open_file($file);
    return $self;
}

sub open_file
{
    my ($self,$file) = @_;

    if( $file ) {
      open(FILE,">$file") or die "Cannot open file: $file\n";
      $self->{file_handle} = \*FILE;
      $self->{file_name} = $file;
    } else {
      $self->{file_handle} = \*STDOUT;
      $self->{file_name} = "-";
    }
}

sub close_file
{
    my ($self) = @_;

    my $file_handle = $self->{file_handle};
    my $file_name   = $self->{file_name};

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

sub write_header {

  my ($self,$title) = @_;

  my $FH = $self->{file_handle};

  print $FH "<html>\n";
  print $FH "<head>\n";
  print $FH "<title>$title</title>\n";
  print $FH "</head>\n";
  print $FH "<body>\n";
  print $FH "\n";
}

sub write_trailer {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "\n";
  print $FH "</body>\n";
  print $FH "</html>\n";
}

sub write_header_redirect {

  my ($self,$title,$html,$seconds) = @_;

# use $html as: url=http://www.yourdomain.com/index.html

  $seconds = 0 unless $seconds;

  my $FH = $self->{file_handle};

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

  my $FH = $self->{file_handle};

  print $FH "\n";
  if( $extra ) {
    print $FH "<table $extra>\n";
  } else {
    print $FH "<table>\n";
  }
}

sub close_table {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "</table>\n";
  print $FH "\n";
}

sub open_table_row {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "<tr>\n";
}

sub close_table_row {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "</tr>\n";
}

sub insert_table_data {

  my ($self,$data,$color) = @_;

  my $FH = $self->{file_handle};

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

  my $FH = $self->{file_handle};

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

  my $FH = $self->{file_handle};

  print $FH "</$type>\n";
  print $FH "\n";
}

sub insert_list {

  my ($self,$text) = @_;

  my $FH = $self->{file_handle};

  print $FH "<li>$text</li>\n";
}

#-------------------------------------------------------------------

sub insert_heading {

  my ($self,$data,$level) = @_;

  my $FH = $self->{file_handle};

  $level = 1 unless $level;
  my $tag = "h" . $level;
  print $FH "<$tag>$data</$tag>\n";
}

sub insert_para {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle};

  print $FH "<p>$data</p>\n";
}

sub insert_data {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle};

  print $FH "$data";
}

sub insert_line {

  my ($self,$data) = @_;

  my $FH = $self->{file_handle};

  print $FH "$data\n";
}

sub insert_break {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "<br>\n";
}

sub insert_image {

  my ($self,$image,$options) = @_;

  my $FH = $self->{file_handle};

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

################################
1;
################################
