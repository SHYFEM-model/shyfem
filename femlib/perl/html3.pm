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
# utilities for reading and writing html files
#
# info:
#
# https://www.w3schools.com/html/default.asp
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use html3;
# use html::tooltip;
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
# version 3.3	01.05.2019	scripts and style elements
# version 3.4	06.06.2020	some enhancements reading tags
# version 3.5	06.06.2020	other enhancements, embed image
# version 3.6	07.09.2020	tooltips and other improvements
# version 3.7	08.10.2020	insert general scripts
# version 3.8	30.11.2020	allow for ancor id
# version 3.9	11.03.2022	avoid parsing empty tag in get_tag()
#
##############################################################

use strict;

package html3;

$::version = "3.8";

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
			 ,version		=>	$::version
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

  my ($self,$title,$rstyle,$extra) = @_;

  $self->write_header_open($title);
  $rstyle->() if $rstyle;
  $self->insert_scripts();
  $self->write_header_close($extra);
}

sub write_header_open {

  my ($self,$title) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<html>\n";
  print $FH "<head>\n";
  print $FH "<title>$title</title>\n";
}

sub write_header_close {

  my ($self,$extra) = @_;

  $extra = "" unless $extra;

  my $FH = $self->{file_handle_write};

  print $FH "</head>\n";
  print $FH "<body $extra>\n";
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

#-------------------------------------------------------------------
# style
#-------------------------------------------------------------------

sub open_style {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "\n";
  print $FH "<style>\n";
}

sub close_style {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "</style>\n";
  print $FH "\n";
}

sub insert_style {

  my ($self,$tag,@values) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "$tag {\n";
  foreach my $value (@values) {
    print $FH "  $value;\n";
  }
  print $FH "}\n";
}

#-------------------------------------------------------------------
# table
#-------------------------------------------------------------------

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

sub insert_table_caption {

  my ($self,$caption) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<caption>$caption</caption>\n";
}

sub open_table_row {

  my ($self,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<tr $extra>\n";
  } else {
    print $FH "<tr>\n";
  }
}

sub close_table_row {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "</tr>\n";
}

sub insert_table_data {

  my ($self,$data,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<td $extra>$data</td>\n";
    #print $FH "<td><font color=\"$color\">$data</font></td>\n";
  } else {
    print $FH "<td>$data</td>\n";
  }
}

sub insert_table_header {

  my ($self,$data,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<th $extra>$data</th>\n";
    #print $FH "<th><font color=\"$color\">$data</font></th>\n";
  } else {
    print $FH "<th>$data</th>\n";
  }
}

#-------------------------------------------------------------------
# list
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
#  insert arbitrary tag
#-------------------------------------------------------------------

sub open_tag {

  my ($self,$tag,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<$tag $extra>\n";
  } else {
    print $FH "<$tag>\n";
  }
}

sub close_tag {

  my ($self,$tag) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "</$tag>\n";
}

#-------------------------------------------------------------------
#  insert text
#-------------------------------------------------------------------

sub insert_heading {

  my ($self,$data,$level,$bookmark) = @_;

  my $FH = $self->{file_handle_write};

  $level = 1 unless $level;
  my $tag = "h" . $level;
  $bookmark = " id=\"$bookmark\"" if $bookmark;
  $bookmark = "" unless $bookmark;

  print $FH "<$tag$bookmark>$data</$tag>\n";
}

sub insert_para {

  my ($self,$data,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<p $extra>$data</p>\n";
  } else {
    print $FH "<p>$data</p>\n";
  }
}

sub insert_pre {

  my ($self,$data,$extra) = @_;

  my $FH = $self->{file_handle_write};

  if( $extra ) {
    print $FH "<pre $extra>$data</pre>\n";
  } else {
    print $FH "<pre>$data</pre>\n";
  }
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

sub insert_vertical_space {

  my ($self,$space) = @_;

  my $FH = $self->{file_handle_write};

  print $FH "<div style=\"position:relative; height:$space\"></div>\n";
}

sub insert_image {

  my ($self,$image,$options) = @_;

  my $FH = $self->{file_handle_write};

  $options = "" unless $options;

  print $FH "<img src=\"$image\" $options>";
}

sub embed_image {

  # before using this, convert your image with "base64 image > image.txt"
  # use format to identify image format: jpg, gif, png, ...
  # use image to give name of image file (image.txt)

  my ($self,$image,$format,$options) = @_;

  my $FH = $self->{file_handle_write};

  $options = "" unless $options;

  print $FH "<img $options src=\"data:image/$format;base64,\n";
  open(IMG,"<$image") || die "cannot open file: $image\n";
  while(<IMG>) {
    print $FH "$_";
  }
  close(IMG);
  print $FH "\">";
}

#-------------------------------------------------------------------
# hyperlink
#-------------------------------------------------------------------

sub make_anchor {

  my ($self,$text,$href,$id) = @_;

  if( $id ) {
    return "<a class=\"$id\" href=\"$href\">$text</a>";
  } else {
    return "<a href=\"$href\">$text</a>";
  }
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
    return "";
  }
}

sub skip_next_tag {

  my ($self,$text,$tag) = @_;

  unless( $tag ) {	#skip next tag whatever
    $tag = $self->find_next_tag($text)
  }

  #print "skipping: $text\n";
  my ($content,$options,$rest) = $self->get_tag($tag,$text);
  #print "skipping rest: $rest\n";

  return $rest;
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

  return unless $text;

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

sub show_lines { # shows limited number of lines

  my ($self,$lines,$n) = @_;

  $n = 10 unless $n;
  my @f = split(/\n/,$lines);
  foreach (@f) {
    print "$_\n";
    $n--;
    last unless $n;
  }
}

sub parse_options {

  my ($self,$options) = @_;

  my ($key,$value);
  my %options = ();

  return \%options unless $options;

  $options =~ s/\n/ /g;

  while ( $options =~ /^\s*([\w-]+)\s*=\s*/ ) {
    #print "key: $options --- $1\n";
    $key=$1;
    $options =~ s/^\s*[\w-]+\s*=\s*//;
    $value="";
    if( $options =~ /^(.*?)\s*([\w-]+)\s*=\s*(.*)$/ ) {
      my ($before,$match,$after) = ($1,$2,$3);
      $value = $1;
      #print "before delete: |$options|\n";
      #print "before match after: |$before|$match|$after|\n";
      $options="$match=$after";
    }
    $options{$key} = $value;
  }
  $options{$key} = $options;
  #print "end: $options \n";
  
  foreach $key (keys %options) {
    $options{$key} =~ s/^\"//;
    $options{$key} =~ s/\"$//;
  }
  return \%options;
}

#-------------------------------------------------------------------
# general scripts
#-------------------------------------------------------------------

sub insert_scripts
{
  my ($self) = @_;

  my $FH = $self->{file_handle_write};

  print $FH <<'EOT';
        <script type="text/javascript">
            function openTab(th)
            {
                window.open(th.name,'_blank');
            }
        </script>
EOT

}

#-------------------------------------------------------------------
# tooltips
#-------------------------------------------------------------------

sub make_tooltip {

  my ($self,$text,$tooltip) = @_;

  my $line = "";
  $line .= "<div class=\"tooltip\">$text\n";
  $line .= "  <span class=\"tooltiptext tooltip-left\">$tooltip</span>\n";
  $line .= "</div>\n";

  return $line;
}

sub init_tooltip {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

print $FH <<'EOT';

<style>
/* Tooltip container */
.tooltip {
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black; /* If you want dots under the hoverable text */
}

/* Tooltip text */
.tooltip .tooltiptext {
  visibility: hidden;
  width: 240px;
  background-color: black;
  color: #fff;
  text-align: center;
  padding: 5px 0;
  border-radius: 6px;
 
  /* Position the tooltip text - see examples below! */
  position: absolute;
  z-index: 1;
}

.tooltip-right {
  top: -5px;
  left: 125%;  
}

.tooltip-left {
  top: -5px;
  bottom:auto;
  right: 128%;  
}

.tooltip-bottom {
  top: 135%;
  left: 50%;  
  margin-left: -60px;
}

.tooltip-top {
  bottom: 125%;
  left: 50%;  
  margin-left: -60px;
}

/* Show the tooltip text when you mouse over the tooltip container */
.tooltip:hover .tooltiptext {
  visibility: visible;
}
</style>

EOT

}

#-------------------------------------------------------------------
# show/hide elements
#-------------------------------------------------------------------

$::pageIdunique = 10000;

sub make_showhide {

  my ($self,$text,$page,$pageId) = @_;

  unless( $pageId ) {
    $pageId = "showhide_id_" . $::pageIdunique++
  }

  #print "id: $pageId\n";
  
  my $line = "";
  $line .= "<button onclick=\"my_showhide(\'$pageId\')\">$text</button>\n";
  $line .= "<div id=\"$pageId\" style=\"display:none\">$page</div>\n";
  #$line .= "<div id=\"$pageId\" style=\"display:block\">$page</div>\n";

  return $line;
}

sub init_showhide {

  my ($self) = @_;

  my $FH = $self->{file_handle_write};

print $FH <<'EOT';

<script>

function my_showhide(pageID) {
  var x = document.getElementById(pageID);
  if (x.style.display === "none") {
    x.style.display = "block";
  } else {
    x.style.display = "none";
  }
} 

</script>

EOT

}

#-------------------------------------------------------------------
# slideshow
#-------------------------------------------------------------------

sub insert_slideshow
{
  my ($self,$imgs,$captions) = @_;

  my $ntot = @$imgs;

  unless( $captions) {
    my @aux = @$imgs;
    $captions = \@aux;
  }

  my $FH = $self->{file_handle_write};

  print $FH "<div class='slideshow-container'>\n";

  for( my $i=1; $i<=$ntot; $i++ ) {
    my $img = shift(@$imgs);
    my $caption = shift(@$captions);
    print $FH "<div class='mySlides fade'>\n";
    #print $FH "  <div class='numbertext'>$i / $ntot</div>\n";
    print $FH "  <img src=\"$img\" style='width:100%'>\n";
    #print $FH "  <div class='text'>$caption</div>\n";
    print $FH "</div>\n";
  }

  print $FH "<a class='prev' onclick='plusSlides(-1)'>&#10094;</a>\n";
  print $FH "<a class='next' onclick='plusSlides(1)'>&#10095;</a>\n";
  print $FH "</div>\n";

  print $FH "<div style=\"text-align:center\">\n";
  #print $FH "<span class=\"dot\" onclick=\"plusSlide(-1)\">-</span>\n";
  for( my $i=1; $i<=$ntot; $i++ ) {
    print $FH "<span class=\"dot\" onclick=\"currentSlide($i)\"></span>\n";
  }
  #print $FH "<span class=\"dot\" onclick=\"plusSlide(1)\">+</span>\n";
  print $FH "</div>\n";
}

sub init_slideshow
{
  my ($self) = @_;

  my $FH = $self->{file_handle_write};

print $FH <<'EOT';

<style>
* {box-sizing: border-box}
body {font-family: Verdana, sans-serif; margin:0}
.mySlides {display: none}
img {vertical-align: middle;}

/* Slideshow container */
.slideshow-container {
  max-width: 1000px;
  position: relative;
  margin: auto;
}

/* Next & previous buttons */
.prev, .next {
  cursor: pointer;
  position: absolute;
  top: 50%;
  width: auto;
  padding: 16px;
  margin-top: -22px;
  color: white;
  font-weight: bold;
  font-size: 18px;
  transition: 0.6s ease;
  border-radius: 0 3px 3px 0;
  user-select: none;
}

/* Position the "next button" to the right */
.next {
  right: 0;
  border-radius: 3px 0 0 3px;
}

/* On hover, add a black background color with a little bit see-through */
.prev:hover, .next:hover {
  background-color: rgba(0,0,0,0.8);
}

/* Caption text */
.text {
  #color: #f2f2f2;
  color: black
  font-size: 15px;
  padding: 8px 12px;
  position: absolute;
  bottom: 8px;
  width: 100%;
  text-align: center;
}

/* Number text (1/3 etc) */
.numbertext {
  #color: #f2f2f2;
  color: black
  font-size: 12px;
  padding: 8px 12px;
  position: absolute;
  top: 0;
}

/* The dots/bullets/indicators */
.dot {
  cursor: pointer;
  height: 15px;
  width: 15px;
  margin: 0 2px;
  background-color: #bbb;
  border-radius: 50%;
  display: inline-block;
  transition: background-color 0.6s ease;
}

.active, .dot:hover {
  background-color: #717171;
}

/* Fading animation */
.fade {
  -webkit-animation-name: fade;
  -webkit-animation-duration: 1.5s;
  animation-name: fade;
  animation-duration: 1.5s;
}

@-webkit-keyframes fade {
  from {opacity: .4} 
  to {opacity: 1}
}

@keyframes fade {
  from {opacity: .4} 
  to {opacity: 1}
}

/* On smaller screens, decrease text size */
@media only screen and (max-width: 300px) {
  .prev, .next,.text {font-size: 11px}
}
</style>

EOT
}

sub init_slideshow_script
{
  my ($self) = @_;

  my $FH = $self->{file_handle_write};

print $FH <<'EOT';

<script>
var slideIndex = 1;
showSlides(slideIndex);

function plusSlides(n) {
  showSlides(slideIndex += n);
}

function currentSlide(n) {
  showSlides(slideIndex = n);
}

function showSlides(n) {
  var i;
  var slides = document.getElementsByClassName("mySlides");
  var dots = document.getElementsByClassName("dot");
  if (n > slides.length) {slideIndex = 1}    
  if (n < 1) {slideIndex = slides.length}
  for (i = 0; i < slides.length; i++) {
      slides[i].style.display = "none";  
  }
  for (i = 0; i < dots.length; i++) {
      dots[i].className = dots[i].className.replace(" active", "");
  }
  slides[slideIndex-1].style.display = "block";  
  dots[slideIndex-1].className += " active";
}
</script>

EOT
}

sub init_meta
{
  my ($self) = @_;

  my $FH = $self->{file_handle_write};

print $FH <<'EOT';
  <meta name="viewport" content="width=device-width, initial-scale=1">
EOT
}

#-------------------------------------------------------------------

################################
1;
################################
