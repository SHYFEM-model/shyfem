#!/usr/bin/perl -w
#
# html utilities
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

###############################################################################

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

sub insert_break {

  my ($self) = @_;

  my $FH = $self->{file_handle};

  print $FH "<br>\n";
}

sub make_anchor {

  my ($self,$text,$href) = @_;

  my $line = "<a href=\"$href\">$text</a>";

  return $line;
}

################################
1;
################################
