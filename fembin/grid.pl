#!/usr/bin/perl -ws
#
# digitize an image (GIF)

use strict;
use Tk;
use Tk::Dialog;
require Tk::FileSelect;

#--------------------------------------------------------
# use grd;			#privat use
# use grdline;			#privat use

# require "image.pl";		#privat use
# require "util_grid.pl";		#privat use
# require "plot.pl";		#privat use
# require "utility.pl";		#privat use
# require "command.pl";		#privat use
# require "open_close.pl";	#privat use
#--------------------------------------------------------

$::image = 1 if $::image and $::image eq "1";

%::appl = ();

my ($cw,$ch) = (800,800);
set_screen(0,0,1);
set_scale(0,0,0,0,1);
tentative_mode("TENT_NODE");
$::appl{canvas_size} = [0,0];

my $image_file = $::image;
my $grid_file = $ARGV[0];

my $mw = MainWindow->new; $mw->title("grid"); $mw->minsize($cw,$ch);
$mw->bind('<Destroy>',\&exec_file_exit);
$::appl{mw} = $mw;

Make_Menu($mw);
Make_Message_bar($mw);

if( $grid_file ) {
  my $error = open_file($grid_file);
  die "Cannot read file : $grid_file\n" if $error;
  $mw->title("grid $grid_file");
} else {
  new_file();
}

read_Image($image_file);
Make_Canvas($mw);
Test_events();

MainLoop();

Write_grid();

#--------------------------------------------------------------------

sub Test_events
{
	my $mw = $::appl{mw};
	#$mw->bind('<Map>',[\&Event,"map"]);
	$mw->bind('<Configure>',[\&Conf,"conf"]);
	#$mw->bind('<ButtonRelease>',[\&Event,"button_release"]);
	#$mw->bind('<Resize>',[\&Event,"resize"]);
}

sub Event
{
	my $count = $::appl{count}++;
	print STDERR "event ($count): @_\n";
}

sub Map
{
	print STDERR "mapping event: @_\n";
}
sub Conf
{
	my $canvas = $::appl{canvas};

	my $count = $::appl{configs}++;
	#print STDERR "Conf configure event ($count): @_\n";
	unless ( $count ) {	#first...
	  print STDERR "Enter first configured...\n";
	  $canvas->Tk::bind('<Enter>',[\&Reconf,"reconf"]);
	}
}

sub Reconf
{
	my $canvas = $::appl{canvas};

	my ($cw_old,$ch_old) = @{$::appl{canvas_size}};
	my ($cw,$ch) = get_Canvas_Size();
	print STDERR "canvas: $cw_old $ch_old  ->  $cw $ch\n";

	$::appl{configs} = 0;
	$canvas->Tk::bind('<Enter>',"");

	if( $cw_old != $cw or $ch_old != $ch ) {
	  print STDERR "re-configure canvas: $cw_old $ch_old  ->  $cw $ch\n";
	  init_file();
	  $::appl{canvas_size} = [$cw,$ch];;
	}
}

#--------------------------------------------------------------------

sub Make_Menu
{
  my ($mw) = @_;

  my $mbar = $mw->Frame(-relief=>'ridge', -borderwidth=>2)
		->pack(-side=>'top', -anchor=>'n', -expand=>0, -fill=>'x');

  my $file = $mbar->Menubutton(-text=>"File", -menuitems => [
		 [ 'command' => "New", -command => \&do_file_new ]
		 ,[ 'command' => "Open", -command => \&do_file_open ]
		 ,[ 'command' => "Close", -command => \&do_file_close ]
		,"-"
		 ,[ 'command' => "Save", -command => \&do_file_save ]
		 ,[ 'command' => "Save as", -command => \&do_file_save_as ]
		,"-"
		,[ 'command' => "Plot", -command => \&do_file_plot ]
		,[ 'command' => "Info", -command => \&do_file_info ]
		,[ 'command' => "Reset", -command => \&do_file_reset ]
		,"-"
		,[ 'command' => "Exit", -command => \&do_file_exit ]
				 ])->pack(-side=>"left");

  my $view = $mbar->Menubutton(-text=>"View", -menuitems => [
		 [ 'command' => "Original", -command => \&do_view_original ]
		,[ 'command' => "Total", -command => \&do_view_total ]
		,"-"
		,[ 'command' => "Zoom in", -command => \&do_view_zoom_in ]
		,[ 'command' => "Zoom out", -command => \&do_view_zoom_out ]
		,[ 'command' => "Window", -command => \&do_view_window ]
		,"-"
		,[ 'command' => "Move", -command => \&do_view_move ]
				 ])->pack(-side=>"left");

  my $mode = $mbar->Menubutton(-text=>"Mode", -menuitems => [
		 [ 'command' => "Node", -command => \&do_mode_node ]
		,[ 'command' => "Element", -command => \&do_mode_elem ]
		,[ 'command' => "Line", -command => \&do_mode_line ]
				 ])->pack(-side=>"left");

  my $node = $mbar->Menubutton(-text=>"Node", -menuitems => [
		 [ 'command' => "Make", -command => \&do_node_make ]
		,[ 'command' => "Delete", -command => \&do_node_delete ]
		,[ 'command' => "Move", -command => \&do_node_move ]
		,[ 'command' => "Unify", -command => \&do_node_unify ]
				 ])->pack(-side=>"left");

  my $elem = $mbar->Menubutton(-text=>"Element", -menuitems => [
		 [ 'command' => "Make", -command => \&do_elem_make ]
		,[ 'command' => "Delete", -command => \&do_elem_delete ]
		,[ 'command' => "Remove", -command => \&do_elem_remove ]
				 ])->pack(-side=>"left");

  my $line = $mbar->Menubutton(-text=>"Line", -menuitems => [
		 [ 'command' => "Make", -command => \&do_line_make ]
		,[ 'command' => "Delete", -command => \&do_line_delete ]
		,[ 'command' => "Remove", -command => \&do_line_remove ]
		,[ 'command' => "Split", -command => \&do_line_split ]
		,[ 'command' => "Join", -command => \&do_line_join ]
		,"-"
		,[ 'command' => "Join all", -command => \&do_line_join_all ]
		,[ 'command' => "Move", -command => \&do_line_move ]
		,[ 'command' => "Close", -command => \&do_line_close ]
				 ])->pack(-side=>"left");

  $::appl{line_end} = 1;
  $::appl{plot_node} = 0;
  $::appl{auto_close} = 0;

  my $options = $mbar->Menubutton(-text=>"Options")->pack(-side=>"left");
		$options->checkbutton(-label=>"Plot line end",
			-command=>\&do_option_line_end,
			-variable=>\$::appl{line_end}
			);
		$options->checkbutton(-label=>"Plot node",
			-command=>\&do_option_plot_node,
			-variable=>\$::appl{plot_node}
			);
		$options->separator();
		$options->checkbutton(-label=>"Autoclose line",
			-command=>\&do_option_auto_close,
			-variable=>\$::appl{auto_close}
			);
		$options->separator();

  Make_Motion_Menu($options);

  my $help = $mbar->Menubutton(-text=>"Help", -menuitems => [
		 [ 'command' => "Help", -command => \&do_nothing ]
		,[ 'command' => "About", -command => \&do_nothing ]
				 ])->pack(-side=>"right");
}

sub Make_Message_bar
{
	my $mw = shift;

	my $message = $mw->Frame()->pack(-side=>'bottom'
			, -anchor=>'s', -expand=>0, -fill=>'x');
	my $mlabel = $message->Label(-text=>"",-relief=>'ridge')
			->pack(-expand=>1, -fill=>'x');

	$::appl{message_bar} = $mlabel;
}

sub Make_Canvas
{
	my $mw = shift;

	my $canvas = $mw->Canvas(-relief=>"sunken", -cursor=>"crosshair", 
		-width=>$cw, -height=>$ch
		)->pack(-side=>'left', -anchor=>'n'
		, -expand=>1, -fill=>'both');
	$canvas->Tk::bind('<Button>'
		,[\&Button_Click,"Click",Ev("b"),Ev("x"),Ev("y")]);

	$::appl{canvas} = $canvas;
}

#-----------------------------------------------------------------

sub do_file_new {	Exec(\&exec_file_new);	}
sub do_file_open {	Exec(\&exec_file_open);	}
sub do_file_close {	Exec(\&exec_file_close);	}
sub do_file_save {	Exec(\&exec_file_save);	}
sub do_file_save_as {	Exec(\&exec_file_save_as);	}
sub do_file_reset {	Dispatch(\&exec_file_reset);	}
sub do_file_plot {	Dispatch(\&exec_file_plot);	}
sub do_file_info {	Dispatch(\&exec_file_info);	}
sub do_file_exit {	Dispatch(\&exec_file_exit);	}

sub do_view_original {	Dispatch(\&exec_view_original);	}
sub do_view_total {	Dispatch(\&exec_view_total);	}
sub do_view_zoom_in {	Dispatch(\&exec_view_zoom_in);	}
sub do_view_zoom_out {	Dispatch(\&exec_view_zoom_out);	}
sub do_view_window {	Dispatch(\&exec_view_window);	}
sub do_view_move {	Dispatch(\&exec_view_move);	}

sub do_mode_node {	Dispatch(\&exec_mode_node);	}
sub do_mode_elem {	Dispatch(\&exec_mode_elem);	}
sub do_mode_line {	Dispatch(\&exec_mode_line);	}

sub do_node_make {	Dispatch(\&exec_node_make);	}
sub do_node_delete {	Dispatch(\&exec_node_delete);	}
sub do_node_move {	Dispatch(\&exec_node_move);	}
sub do_node_unify {	Dispatch(\&exec_node_unify);	}

sub do_elem_make {	Dispatch(\&exec_elem_make);	}
sub do_elem_delete {	Dispatch(\&exec_elem_delete);	}
sub do_elem_remove {	Dispatch(\&exec_elem_remove);	}

sub do_line_make {	Dispatch(\&exec_line_make);	}
sub do_line_delete {	Dispatch(\&exec_line_delete);	}
sub do_line_remove {	Dispatch(\&exec_line_remove);	}
sub do_line_split {	Dispatch(\&exec_line_split);	}
sub do_line_join {	Dispatch(\&exec_line_join);	}
sub do_line_join_all {	Dispatch(\&exec_line_join_all);	}
sub do_line_move {	Dispatch(\&exec_not_yet);	}
sub do_line_close {	Dispatch(\&exec_line_close);	}

sub do_option_line_end {	Dispatch(\&exec_option_line_end);	}
sub do_option_plot_node {	Dispatch(\&exec_option_plot_node);	}
sub do_option_auto_close {	Exec(\&exec_option_auto_close);	}
sub do_option_motion {		Exec(\&exec_option_motion);	}

sub do_nothing {	Dispatch(\&exec_nothing);	}

#-----------------------------------------------------------------

sub act_command
{
	$::appl{act_command} = $_[0] if defined $_[0];
	return $::appl{act_command};
}

sub Exec
{
	my $command = $_[0];

	&$command("start") if $command;
}

sub Dispatch
{
	my $command = act_command();

	&$command("end") if $command;
	$command = act_command($_[0]);
	&$command("start") if $command;
}

sub Button_Click
{
	my ($w,$t,$b,$x,$y) = @_;

	my $command = act_command();

	if( $b == 1 ) {
	  &$command($w,$t,$b,$x,$y) if $command;
	} elsif( $b == 3 ) {
	  get_tentative($x,$y);
	} else {
	  print STDERR "*** button $b not linked...\n";
	}
}

sub Bell
{
	my $canvas = $::appl{canvas};
	$canvas->bell();
	print STDERR "Bell called\n";
}

sub Make_Motion_Menu
{
	my $menub = shift;

	$::appl{motion_speed} = 0;
	$::appl{motion_events} = 0;
	$::appl{motion_on} = 0;
	$::appl{motion_command} = "";

	my $motion = $menub->menu->Menu();

	foreach my $val (qw/0 1 2 3 5 8 10 20 30 40 50 80 100/) {
	  $motion->radiobutton(-label=>$val, 
			-command => \&do_option_motion,
	                -variable => \$::appl{motion_speed}, 
			-value => $val );
	}

	$menub->cascade(-label=>"Motion");
	$menub->entryconfigure("Motion",-menu=>$motion);
}

#------------------------------------------------------------

sub tentative_mode
{
	$::appl{act_tentative_mode} = $_[0] if defined $_[0];
	return $::appl{act_tentative_mode};
}

sub reset_tentative_mode
{
	reset_tentative_item();
	$::appl{act_tentative_mode} = $_[0] if defined $_[0];
	return $::appl{act_tentative_mode};
}

sub tentative_item
{
	$::appl{act_tentative_item} = $_[0] if defined $_[0];
	return $::appl{act_tentative_item};
}

sub tentative_extra
{
	my @f = @_;

	$::appl{act_tentative_extra} = \@f if defined $f[0];
	return $::appl{act_tentative_extra};
}

sub reset_tentative_item
{
	tentative_item("");
}

sub info_tentative
{
	my $text = shift;
	$text = "" unless defined $text;

        my $mode = tentative_mode();
        my $item = tentative_item();

	print STDERR "info_tentative ($text): $mode  $item\n";
}

sub write_message_bar
{
	my $text = shift;

	my $message = $::appl{message_bar};
	$message->configure(-text=>"$text");

	print STDERR "message bar: $text\n";
}

sub get_tentative
{
	my ($x,$y) = @_;

	my $maxdist = 100;
	my $mode = tentative_mode();
        my $item = tentative_item();
        #my $grid = $::appl{grid};

	my ($xr,$yr) = S2R($x,$y);

	print STDERR "tentative point $x $y  ->  $xr $yr\n";

	if( $mode eq "TENT_NODE" ) {
	  my $item = get_closest_node($xr,$yr,$item);
	  $item = check_max_dist($item,$x,$y);
	  tentative_item($item);
	  if( $item ) {
	    write_message_bar("Node   $item->{number}    $item->{type}    " . 
			"$item->{x}    $item->{y}");
	  }
	} elsif( $mode =~ "TENT_LINE" ) {
	  my ($item,$xc,$yc,$node,$seg) = get_closest_line($xr,$yr,$item);
	  $item = check_max_dist($item,$x,$y,$xc,$yc);
	  print STDERR "setting tentative line: $mode\n";
	  tentative_item($item);
	  tentative_extra($xc,$yc,$node,$seg);
	  info_tentative("get_tentative: line");
	  if( $item ) {
	    write_message_bar("Line   $item->{number}    $item->{type}    " . 
			"$item->{nvert}");
	  }
	} elsif( $mode =~ "TENT_ELEM" ) {
	  my ($item,$xc,$yc,$node,$seg) = get_closest_elem($xr,$yr,$item);
	  print STDERR "setting tentative elem: $mode\n";
	  #$grid->elem_info($item,"tent elem: ");
	  tentative_item($item);
	  tentative_extra($xc,$yc,$node,$seg);
	  info_tentative("get_tentative: elem");
	  if( $item ) {
	    write_message_bar("Element   $item->{number}    $item->{type}    " . 
			"$item->{nvert}");
	  }
	} else {
	  print STDERR "*** unknown tentative mode: $mode\n";
	}

	write_message_bar("") unless tentative_item();

	Plot_tentative();
}

sub check_max_dist
{
	my ($item,$x,$y,$xc,$yc) = @_;

	my $maxdist = 100;

	if( $item ) {
	    unless( $xc ) {
	      ($xc,$yc) = ($item->{x},$item->{y});
	    }
	    my ($xs,$ys) = R2S($xc,$yc);
	    my $dist = dist($x,$y,$xs,$ys);
	    $item = "" if( $dist > $maxdist );
	}

	return $item;
}

sub make_act_node
{
	my ($b,$x,$y) = @_;

        my $item = get_active_node();

        unless( $item ) {
          my ($xr,$yr) = S2R($x,$y);
          $item = make_node($b,$xr,$yr);
        }
        Plot_node($item,-tags=>"ActItem");

	return $item;
}

sub motion_speed
{
	return $::appl{motion_speed};
}

#----------------------------------------

sub get_global
{
	return $::appl{$_[0]};
}

sub set_global
{
	my ($what,$value) = @_;

	$::appl{$what} = $value;
}

sub set_and_get_global
{
	my ($what,$value) = @_;

	my $ref = \$::appl{$what};
        $$ref = $value if defined $value;
        return $$ref;
}

sub set_grid_and_file
{
	my ($grid,$file) = @_;

        $::appl{grid} = $grid;
        $::appl{act_grid_file} = $file;
}

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------

sub get_active_node
{
        my $mode = tentative_mode();
        my $item = tentative_item();

	if( $mode eq "TENT_NODE" ) {
	  return $item;
	}
}

sub get_active_elem
{
        my $mode = tentative_mode();
        my $item = tentative_item();

	if( $mode =~ "TENT_ELEM" ) {
	  return $item;
	}
}

sub get_active_line
{
        my $mode = tentative_mode();
        my $item = tentative_item();

	if( $mode =~ "TENT_LINE" ) {
	  return $item;
	}
}

#-----------------------------------------------------------------

sub clean_item_list
{
	$::appl{node_list} = [];
}

sub add_to_item_list
{
	my ($item) = @_;

	push(@{$::appl{node_list}},$item);
}

sub fill_of_item_list
{
	return scalar(@{$::appl{node_list}});
}

sub get_item_list
{
	return $::appl{node_list};
}

sub get_item_of_item_list
{
	my ($i) = @_;

	my $list = $::appl{node_list};

	return $list->[$i];
}

sub get_first_of_item_list
{
	my $list = $::appl{node_list};
	my $n = @$list;

	return if $n < 1;
	return $list->[0];
}

sub get_last_of_item_list
{
	my $list = $::appl{node_list};
	my $n = @$list;

	return if $n < 1;
	return $list->[$n-1];
}

sub identical_item_list
{
	my $list = $::appl{node_list};
	my $n = @$list;

	return 0 if $n < 2;
	return 0 if $list->[$n-1] ne $list->[$n-2];

	return 1;
}

#-----------------------------------------------------------------

sub Set_motion
{
        my $command = shift;

        my $canvas = $::appl{canvas};

        if( $::appl{motion_on} ) {
            $canvas->Tk::bind('<Motion>',"");
        } else {
            $::appl{motion_events} = 0;
            $canvas->Tk::bind('<Motion>'
                        ,[\&Handle_motion,"Motion",0,Ev("x"),Ev("y")]);
        }

        $::appl{motion_on} = toggle( $::appl{motion_on} );
        $::appl{motion_command} = $command;

        return $::appl{motion_on};
}

sub Handle_motion
{
        my ($w,$text,$b,$x,$y) = @_;

        my $motion_speed = $::appl{motion_speed};
        my $motion_events = $::appl{motion_events};
        my $motion_command = $::appl{motion_command};

        if( $motion_speed <= 0 ) {
          print STDERR "*** error in Handle_motion: $motion_speed\n";
        } elsif( $motion_events % $motion_speed == 0 ) {
          &$motion_command($w,$text,$b,$x,$y);
        }

        $::appl{motion_events}++;
}

sub Handle_drag
{
        my ($w,$t,$b,$x,$y,$extra) = @_;

        my $canvas = $::appl{canvas};

        if( $t eq "Click" ) {
            print STDERR "Handle_drag - button click: $t  $b  $x  $y\n";
            $canvas->Tk::bind('<ButtonRelease>',[\&Handle_drag,"Release"
                                ,0,Ev("x"),Ev("y")]);
            $canvas->Tk::bind('<Motion>',[\&Handle_drag,"Motion"
                                ,0,Ev("x"),Ev("y")]);

            $::appl{drag_xy} = [$x,$y];
            $::appl{drag_command} = $extra;
        } elsif( $t eq "Release" ) {
            print STDERR "Handle_drag - button release: $t  $b  $x  $y\n";
            $canvas->Tk::bind('<ButtonRelease>',"");
            $canvas->Tk::bind('<Motion>',"");
            Delete_act();
            my $xyold = $::appl{drag_xy};
	    my $command = $::appl{drag_command};
	    &$command("area","area",$b,$x,$y,$xyold);
        } elsif( $t eq "Motion" ) {
            print STDERR "Handle_drag - motion: $t  $b  $x  $y\n";
            my ($xold,$yold) = @{$::appl{drag_xy}};
	    $canvas->delete("ActItem");
            $canvas->createRectangle($xold,$yold,$x,$y,-tags=>"ActItem");
        } else {
            print STDERR "Handle_drag - impossible: $t  $b  $x  $y\n";
        }
}

#-----------------------------------------------------------------

sub choose_file
{
	my ($defext,$file) = @_;

        my $mw = $::appl{mw};
	my $types = [ ['My Files',       ["$defext"]] ];

	if( defined $file ) {
	  $file = $mw->getSaveFile( -initialfile => $file 
			, -filetypes => $types );
	} else {
	  $file = $mw->getOpenFile( -filetypes=>$types );
	}

	$file = strip_dir($file);
	return $file;
}

sub strip_dir
{
	my $path = shift;

	$path =~ s/^.*\/// if $path;

	return $path;
}

#-----------------------------------------------------------------



sub Dialog_Info_ok
{
	my ($text,$options) = @_;

	my $mw = $::appl{mw};

	my $dialog = $mw->Dialog(-text => "$text", 
		-bitmap => 'info', -title => 'Information', 
		-default_button => 'Ok', -buttons => [qw/Ok/]);

	#my $answer = $dialog->Show($options);
	my $answer = $dialog->Show(-global);
}

sub Dialog_YesNoCancel
{
	my ($text,$options) = @_;

	my $mw = $::appl{mw};

	my $dialog = $mw->Dialog(-text => "$text", 
		-bitmap => 'question', -title => 'Question', 
		-default_button => 'Yes', -buttons => [qw/Yes No Cancel/]);

	return $dialog->Show(-global);
}

# options may be "-global"
#
#	$dialog = $mw->Dialog(-text => 'Save File?', 
#		-bitmap => 'question', -title => 'Save File Dialog', 
#		-default_button => 'Yes', -buttons => [qw/Yes No Cancel/]);
#
# bitmaps: 	(see also /usr/lib/tk8.4/msgbox.tcl)
#
#	error info question default
#
#        switch -- $data(-icon) {
#            "error"     {set data(-icon) "stop"}
#            "warning"   {set data(-icon) "caution"}
#            "info"      {set data(-icon) "note"}
#        }
#
#    switch -- $data(-type) {
#        abortretryignore {
#            set names [list abort retry ignore]
#            set labels [list &Abort &Retry &Ignore]
#        }
#        ok {
#            set names [list ok]
#            set labels {&OK}
#        }
#        okcancel {
#            set names [list ok cancel]
#            set labels [list &OK &Cancel]
#        }
#        retrycancel {
#            set names [list retry cancel]
#            set labels [list &Retry &Cancel]
#        }
#        yesno {
#            set names [list yes no]
#            set labels [list &Yes &No]
#        }
#        yesnocancel {
#            set names [list yes no cancel]
#            set labels [list &Yes &No &Cancel]
#        }
#        default {
#            error "bad -type value \"$data(-type)\": must be\
#                    abortretryignore, ok, okcancel, retrycancel,\
#                    yesno, or yesnocancel"
#
#-----------------------------------------------------------------

#!/usr/bin/perl -w

use strict;

#-----------------------------------------------------------------

sub exec_nothing
{
	my ($w,$t,$b,$x,$y) = @_;
	print STDERR "Nothing: $w\n";

	Dialog_Info_ok("No action");
}

sub exec_not_yet
{
	my ($w,$t,$b,$x,$y) = @_;
	print STDERR "Not yet: $w\n";

	Dialog_Info_ok("Not yet implemented");
}

sub exec_file_new
{
	New_file();
}

sub exec_file_open
{
	Open_file();
}

sub exec_file_close
{
	Close_file();
}

sub exec_file_save
{
	Save_file();
}

sub exec_file_save_as
{
	Save_as_file();
}

sub exec_file_reset
{
	Plot_all();
}

sub exec_file_info
{
	print STDERR "Info\n";
	my ($cw,$ch) = get_Canvas_Size();
	print STDERR "Canvas size is: $cw $ch\n";

	Show_tag_list();
}

sub exec_file_plot
{
	Plot_Postscript();
}

sub exec_file_exit
{
	Exit_file();
}

#-----------------------------------------------------------------

sub exec_view_total
{
	print STDERR "Total View\n";

	my ($pw,$ph) = get_Viewport_Size();
	my ($cw,$ch) = get_Canvas_Size();

	print STDERR "exec_view_total: canvas size is $cw,$ch\n";

	my $zoom = min($cw/$pw,$ch/$ph);
	$zoom = round_binary($zoom) if has_image();
	print STDERR "zoom used: $zoom\n";

	set_new_Viewport();

        my $xc = $pw/2;
        my $yc = $ph/2;
	set_screen(0,0,1);	#we need this to make conversion S2P

	move_view($xc,$yc,$zoom);
	Plot_all();
}

sub exec_view_original
{
	print STDERR "Original View\n";

	my ($pw,$ph) = get_Viewport_Size();

	set_original_Viewport();

        my $xc = $pw/2;
        my $yc = $ph/2;
	my $zoom = 1;
	set_screen(0,0,$zoom);

	move_view($xc,$yc,$zoom);
	Plot_all();
}

sub exec_view_zoom_out
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Zoom out: $w\n";
	  return;
	}

	print STDERR "Zoom out: $t  $b  $x  $y\n";

	my $zoom = get_screen_zoom() / 2;
	move_view($x,$y,$zoom);
	Plot_all();
}

sub exec_view_zoom_in
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Zoom in: $w\n";
	  return;
	}

	print STDERR "Zoom in: $t  $b  $x  $y\n";

	my $zoom = get_screen_zoom() * 2;
	move_view($x,$y,$zoom);
	Plot_all();
}

sub exec_view_move
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "View move: $w\n";
	  return;
	}

	print STDERR "View move: $t  $b  $x  $y\n";

	my $zoom = get_screen_zoom();
	move_view($x,$y,$zoom);
	Plot_all();
}

sub exec_view_window
{
	my ($w,$t,$b,$x,$y,$extra) = @_;

	if( not ref($w) ) {
	  if( $w eq "start" ) {
	    print STDERR "View window: $w\n";
	    Delete_act();
	  } elsif( $w eq "end" ) {
	    print STDERR "View window: $w\n";
	    Delete_act();
	  } elsif( $w eq "area" ) {
	    print STDERR "View window: $w\n";
	    Delete_act();
	    my ($xold,$yold) = @$extra;
	    my ($xc,$yc) = (($x+$xold)/2,($y+$yold)/2);
	    my $zoom = get_screen_zoom();
	    my $dx = abs($x-$xold)/$zoom;
	    my $dy = abs($y-$yold)/$zoom;
	    my ($cw,$ch) = get_Canvas_Size();
	    $zoom = min($cw/$dx,$ch/$dy);
	    $zoom = round_binary($zoom) if has_image();
	    move_view($xc,$yc,$zoom);
	    Plot_all();
	    print STDERR "zoom used: $zoom\n";
	  } else {
	    print STDERR "View window: impossible   $w\n";
	  }
	} else {

	  Handle_drag($w,$t,$b,$x,$y,\&exec_view_window);

	}
}

#--------------------------------------------------------------------

sub exec_mode_node
{
	my ($w) = @_;
	if( not ref($w) and $w eq "start" ) {
	  reset_tentative_mode("TENT_NODE");
	}
}

sub exec_mode_elem
{
	my ($w) = @_;
	if( not ref($w) and $w eq "start" ) {
	  reset_tentative_mode("TENT_ELEM");
	}
}

sub exec_mode_line
{
	my ($w) = @_;
	if( not ref($w) and $w eq "start" ) {
	  reset_tentative_mode("TENT_LINE");
	}
}

#--------------------------------------------------------------------

sub exec_node_make
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Node make: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_NODE");
	    Delete_tentative();
	  }
	  return;
	}

	if( motion_speed() and $b ) {	#real mouse click
	  Set_motion(\&exec_node_make);
	  return;
	}

	print STDERR "Node make: $t  $b  $x  $y\n";

	my ($xr,$yr) = S2R($x,$y);

	my $item = get_active_node();
	if ( $item ) {
	  ($xr,$yr) = ($item->{x},$item->{y});
	}
	my $nitem = make_node($b,$xr,$yr);

	Plot_node($nitem);
	Delete_tentative();
}

sub exec_node_delete
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Node delete: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_NODE");
	  }
	  return;
	}

	print STDERR "Node delete: $t  $b  $x  $y\n";

	my $item = get_active_node();
	return unless $item;

	if( is_used($item) ) {
	  print STDERR "*** Cannot delete used node...\n";
	  return;
	}

	Delete_tentative();
	Unplot_node($item);
	delete_node($item);

	print STDERR "Node deleted...\n";
}

sub exec_node_move
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Node move: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_NODE");
	  }
	  return;
	}

	print STDERR "Node move: $t  $b  $x  $y\n";

	my $item = get_active_node();
	return unless $item;

	($item->{x},$item->{y}) = S2R($x,$y);

	Unplot_tentative();
	Replot_node($item);
	Plot_tentative();
}

sub exec_node_unify
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  if( $w eq "start" ) {
	    print STDERR "Node unify: $w\n";
	    clean_item_list();
	    reset_tentative_mode("TENT_NODE");
	  } else {
	    print STDERR "Node unify: impossible   $w\n";
	  }
	  return;
	}

	print STDERR "Node unify: $t  $b  $x  $y\n";

	my $item = get_active_node();
	return unless $item;

	my $n = fill_of_item_list();

	if( $n == 0 ) {
	  add_to_item_list($item);
	  return;
	} elsif( $n > 1 ) {
	  print STDERR "Node unify: error in item list: $n\n";
	}

	my $first_item = get_item_of_item_list(0);	#this one is deleted
	clean_item_list();

	Delete_tentative();
	Unplot_node($first_item);

	unify_nodes($item,$first_item);

	Replot_node($item);
}

#--------------------------------------------------------------------

sub exec_line_make
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line make: $w\n";
	  if( $w eq "start" ) {
	    clean_item_list();
	    reset_tentative_mode("TENT_NODE");
	  } elsif( $w eq "end" ) {
	    my $item = make_line(get_item_list());
	    auto_close_line($item);
	    clean_item_list();
	    Delete_act();
	    Plot_line($item);
	  }
	  return;
	}

	if( motion_speed() and $b ) {	#real mouse click
	  my $in_motion = Set_motion(\&exec_line_make);
	  exec_line_make("end") unless $in_motion;	#finish line
	  return;
	}

	print STDERR "Line make: $t  $b  $x  $y\n";

	my $item = make_act_node($b,$x,$y);
	Delete_tentative();

	my $last_item = get_last_of_item_list();
	if( $last_item and $last_item eq $item ) {	#finish line
	  print STDERR "Line make: finishing line\n";
	  $item = make_line(get_item_list());
	  auto_close_line($item);
	  clean_item_list();
	  Delete_act();
	  Plot_line($item);
	} else {
	  add_to_item_list($item);
	  if( $last_item ) {
	    my ($x1,$y1) = R2S($last_item->{x},$last_item->{y});
	    my ($x2,$y2) = R2S($item->{x},$item->{y});
	    Plot_line_segment($x1,$y1,$x2,$y2,-tags=>"ActItem");
	  }
	}
}

sub exec_line_delete
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line delete: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_LINE");
	  }
	  return;
	}

	print STDERR "Line delete: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

	my @nodes = @{$item->{vert}};

        Delete_tentative();
        Unplot_line($item);
        delete_line($item);
	Plot_single_nodes(@nodes);

        print STDERR "Line deleted...\n";
}

sub exec_line_remove
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line remove: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_LINE");
	  }
	  return;
	}

	print STDERR "Line remove: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

	my @nodes = @{$item->{vert}};

        Delete_tentative();
        Unplot_line($item);
        delete_line($item);
	Unplot_node_list(@nodes);
	delete_single_nodes(@nodes);

        print STDERR "Line removed...\n";
}

sub exec_line_split
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line split: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_LINE_NODE");
	  }
	  return;
	}

	print STDERR "Line split: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

	my $extra = tentative_extra();
	my $node = $extra->[2];

        Delete_tentative();
        Unplot_line($item);
	my ($line1,$line2) = split_line($item,$node);
	Plot_line($line1);
	Plot_line($line2);
}

sub exec_line_join
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line join: $w\n";
	  if( $w eq "start" ) {
	    clean_item_list();
	    reset_tentative_mode("TENT_LINE_NODE");
	  }
	  return;
	}

	print STDERR "Line join: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

	my $n = fill_of_item_list();

	if( $n == 0 ) {
	  add_to_item_list($item);
	  return;
	} elsif( $n > 1 ) {
	  print STDERR "Line join: error in item list: $n\n";
	}

	my $first_item = get_item_of_item_list(0);	#this is first line
	my $extra = tentative_extra();
	my $node = $extra->[2];

	clean_item_list();

	Delete_tentative();
	Unplot_line($first_item);
	Unplot_line($item);

	my ($new_line) = join_line($first_item,$item,$node);

	if( $new_line ) {
	  Plot_line($new_line);
	} else {
	  Plot_line($first_item);
	  Plot_line($item);
	}
}

sub exec_line_join_all
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line join all: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_LINE");
	  }
	  return;
	}

	print STDERR "Line join all: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

        Delete_tentative();
        Unplot_line($item);
	print STDERR "starting join all...\n";
	my $lines = join_all($item);
	Bell() unless scalar(@$lines);
	Unplot_lines($lines);
	delete_lines($lines);
	Plot_line($item);
	print STDERR "end of join_all\n";
}

sub exec_line_close
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Line close: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_LINE");
	  }
	  return;
	}

	print STDERR "Line remove: $t  $b  $x  $y\n";

        my $item = get_active_line();
        return unless $item;

        Delete_tentative();
        Unplot_line($item);
	close_line($item);
        Plot_line($item);

        print STDERR "Line closed...\n";
}

#--------------------------------------------------------------------

sub exec_elem_make
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Element make: $w\n";
	  if( $w eq "start" ) {
	    clean_item_list();
	    reset_tentative_mode("TENT_NODE");
	  } elsif( $w eq "end" ) {
	    my $item = make_elem(get_item_list());
	    clean_item_list();
	    Delete_act();
	    Plot_elem($item);
	  }
	  return;
	}

	if( motion_speed() and $b ) {	#real mouse click
	  my $in_motion = Set_motion(\&exec_elem_make);
	  exec_elem_make("end") unless $in_motion;	#finish elem
	  return;
	}

	print STDERR "Element make: $t  $b  $x  $y\n";

	my $item = make_act_node($b,$x,$y);
	Delete_tentative();

	my $first_item = get_first_of_item_list();
	if( $first_item and $first_item eq $item ) {	#finish elem
	  print STDERR "Element make: finishing elem\n";
	  $item = make_elem(get_item_list());
	  clean_item_list();
	  Delete_act();
	  Plot_elem($item);
	} else {
	  my $last_item = get_last_of_item_list();
	  add_to_item_list($item);
	  if( $last_item ) {
	    my ($x1,$y1) = R2S($last_item->{x},$last_item->{y});
	    my ($x2,$y2) = R2S($item->{x},$item->{y});
	    Plot_line_segment($x1,$y1,$x2,$y2,-tags=>"ActItem");
	  }
	}
}

sub exec_elem_delete
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Element delete: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_ELEM");
	  }
	  return;
	}

	print STDERR "Element delete: $t  $b  $x  $y\n";

        my $item = get_active_elem();
        return unless $item;

	my @nodes = @{$item->{vert}};

        Delete_tentative();
        Unplot_elem($item);
        delete_elem($item);
	Plot_single_nodes(@nodes);

        print STDERR "Element deleted...\n";
}

sub exec_elem_remove
{
	my ($w,$t,$b,$x,$y) = @_;

	if( not ref($w) ) {
	  print STDERR "Element remove: $w\n";
	  if( $w eq "start" ) {
	    reset_tentative_mode("TENT_ELEM");
	  }
	  return;
	}

	print STDERR "Element remove: $t  $b  $x  $y\n";

        my $item = get_active_elem();
        return unless $item;

	my @nodes = @{$item->{vert}};

        Delete_tentative();
        Unplot_elem($item);
        delete_elem($item);
	Unplot_node_list(@nodes);
	delete_single_nodes(@nodes);

        print STDERR "Element removed...\n";
}

#--------------------------------------------------------------------

sub exec_option_line_end
{
	my ($w) = @_;

	return if $w eq "end";

	my $val = get_global("line_end");
	print STDERR "Line end called ($w): $val\n";

	Plot_all();
}

sub exec_option_plot_node
{
	my ($w) = @_;

	return if $w eq "end";

	my $val = get_global("plot_node");
	print STDERR "Plot node called ($w): $val\n";

	Plot_all();
}

sub exec_option_auto_close
{
	my ($w) = @_;

	return if $w eq "end";

	my $val = get_global("auto_close");
	print STDERR "Auto close called ($w): $val\n";
}

sub exec_option_motion
{
	my ($w) = @_;

	return if $w eq "end";

	my $val = get_global("motion_speed");
	if( get_global("motion_on") and $val <= 0 ) {
	  Set_motion();
	}
	print STDERR "motion called ($w): $val\n";
}

#--------------------------------------------------------------------
1;
#--------------------------------------------------------------------

#!/usr/bin/perl

use strict;

sub read_Image
{
	my $image_file = shift;

	my $mw = $::appl{mw};
	my $pic;

	if( $image_file ) {
	  $pic = $mw->Photo( -file => $image_file );
	  my $ph = $pic->height(); 
	  my $pw = $pic->width();
	  print STDERR "readImage: Size of picture: $pw $ph\n";
	}

	$::appl{pic} = $pic;
}

sub make_Image
{
	my $mw = $::appl{mw};
	my $canvas = $::appl{canvas};

	my ($cw,$ch) = get_Canvas_Size();
	my ($pw,$ph) = get_Image_Size();

	print STDERR "Making Image: $cw,$ch,$pw,$ph\n";

	my $xc = int($pw/2);
	my $yc = int($ph/2);
	set_screen(0,0,1);

	my $showpic = $mw->Photo( -width => $cw , -height => $ch );
	$::appl{showpic} = $showpic;

	move_Image($xc,$yc,1);
}

sub move_Image
{
	my ($x,$y,$zoom) = @_;

	print STDERR "move_Image: x,y,zoom = $x $y $zoom\n";

        reposition_Image($x,$y,$zoom);
}

sub get_Viewport_Size 
{
	if( has_image() ) {
	  return get_Image_Size();
	} else {
	  return get_Canvas_Size();
	}
}

sub get_Image_Size
{
	my $pic = $::appl{pic};

	my ($pw,$ph) = ($pic->width(),$pic->height()); 

	return ($pw,$ph);
}

sub get_Canvas_Size
{
	my $canvas = $::appl{canvas};

	#my ($cw,$ch) = ($canvas->cget(-width),$canvas->cget(-height));
	my ($cw,$ch) = ($canvas->width(),$canvas->height());

	return ($cw,$ch);
}

sub reposition_Image
{
	my ($x,$y,$zoom) = @_;

	my $showpic = $::appl{showpic};
	my $pic = $::appl{pic};

	print STDERR "reposition_Image: zoom = $zoom\n";

	my ($xn0,$yn0,$xn1,$yn1) = set_new_screen_coords($x,$y,$zoom);

	print STDERR "reposition_Image: $zoom\n";
	print STDERR "reposition_Image: $xn0,$yn0,$xn1,$yn1\n";

	$showpic->blank();
	return unless $pic;

	if( $zoom > 1 ) {
	  $showpic->copy($pic,-from => $xn0,$yn0,$xn1,$yn1, -zoom=>$zoom);
	} elsif( $zoom < 1 ) {
	  my $r = int(0.5+1/$zoom);
	  $showpic->copy($pic,-from => $xn0,$yn0,$xn1,$yn1, -subsample=>$r);
	} else {
	  $showpic->copy($pic,-from => $xn0,$yn0,$xn1,$yn1);
	}
}

sub set_new_screen_coords
{
	my ($x,$y,$zoom) = @_;

	my ($pw,$ph) = get_Viewport_Size();
	my ($cw,$ch) = get_Canvas_Size();

	my ($xc,$yc) = S2P($x,$y);
	my ($xn0,$yn0,$xn1,$yn1) = new_coords($xc,$yc,$cw,$ch,$pw,$ph,$zoom);

	set_screen($xn0,$yn0,$zoom);

	return ($xn0,$yn0,$xn1,$yn1);
}

sub new_coords		#make (x,y) the center of the image
{
	my ($x,$y,$cx,$cy,$cxmax,$cymax,$zoom) = @_;

	my ($x0,$y0,$x1,$y1);

	if( has_image() ) {
	  ($x0,$x1) = adjust_coord($x,$cx,$cxmax,$zoom);
	  ($y0,$y1) = adjust_coord($y,$cy,$cymax,$zoom);
	} else {
	  ($x0,$x1) = compute_coord($x,$cx,$cxmax,$zoom);
	  ($y0,$y1) = compute_coord($y,$cy,$cymax,$zoom);
	}

	print STDERR "new_coords: $x0,$y0,$x1,$y1 - $cxmax,$cymax\n";

	return ($x0,$y0,$x1,$y1);
}

sub compute_coord
{
	my ($xy,$w,$m,$zoom) = @_;

	$w = $w / $zoom;
	my $xy0 = $xy - $w/2;
	my $xy1 = $xy + $w/2;

	return ($xy0,$xy1);
}

sub adjust_coord
{
	my ($xy,$w,$m,$zoom) = @_;

	my ($xy0,$xy1);
	$w = $w / $zoom;

	if( $w >= $m ) {
	  $xy0 = 0;
	  $xy1 = $m;
	} else {
	  $xy0 = $xy - $w/2;
	  $xy1 = $xy + $w/2;
	  if( $xy0 < 0 ) {
	    $xy0 = 0;
	    $xy1 = $w;
	  } elsif( $xy1 >= $m ) {
	    $xy1 = $m;
	    $xy0 = $xy1-$w;
	  }
	}

	return ($xy0,$xy1);
}

#------------------------------------------------------

sub has_image
{
	return $::appl{pic};
}

sub set_screen
{
	my ($xs0,$ys0,$szoom) = @_;
	$::appl{screen_origin} = [$xs0,$ys0];
	$::appl{szoom} = $szoom;
}

sub set_screen_origin
{
	my ($xs0,$ys0) = @_;
	$::appl{screen_origin} = [$xs0,$ys0];
}

sub get_screen_origin
{
	return @{$::appl{screen_origin}};
}

sub set_screen_zoom
{
	my ($szoom) = @_;
	$::appl{szoom} = $szoom;
}

sub get_screen_zoom
{
	return $::appl{szoom};
}

#------------------------------------------------------

sub S2P
{
	my ($xs,$ys) = @_;

	my ($xs0,$ys0) = @{$::appl{screen_origin}};
	my $szoom = $::appl{szoom};

	my $xp = $xs0 + $xs / $szoom;
	my $yp = $ys0 + $ys / $szoom;

	#return (int $xp,int $yp);
	return ($xp,$yp);
}

sub P2S
{
	my ($xp,$yp) = @_;

	my ($xs0,$ys0) = @{$::appl{screen_origin}};
	my $szoom = $::appl{szoom};

	my $xs = ($xp - $xs0) * $szoom;
	my $ys = ($yp - $ys0) * $szoom;

	#return (int $xs,int $ys);
	return ($xs,$ys);
}

#--------------------------------------
1;
#--------------------------------------

#!/usr/bin/perl -w

use strict;

#------------------------------------------------------------------
#
#
#------------------------------------------------------------------

sub Close_file
{
	if( not is_open() ) {
		return;
	}

	if( is_changed() ) {
		my $answer = Do_you_want_to_save();
		if( $answer eq "Cancel" ) {
			my $error = "File not closed.";
			return $error;
		} elsif( $answer eq "Yes" ) {
			my $error = Save_as_file();
			return $error if $error;
		}
	}

	my $file = act_file();
	my $error = close_file($file);
	return $error if $error;

	reset_file();

	return;
}

sub New_file
{
	if( is_open() ) {
		my $error = Close_file();
		return $error if $error;
	}

	new_file();
	init_file();

	return;
}

sub Open_file
{
	if( is_open() ) {
		my $error = Close_file();
		return $error if $error;
	}

	my $file = get_name();
	if( $file ) {
		my $error = open_file($file);
		return $error if $error;
		init_file($file);
	} else {
		my $error = "No file given.";
		return $error;
	}

	is_changed("");
	is_saved("");

	return;
}

sub Exit_file
{
	my $error = Close_file();

	return $error if $error;

	exit_file();
}

sub Save_file
{
	if( not is_open() or not is_changed() ) {
		return;
	}

	if( is_saved() ) {
		my $file = act_file();
		my $error = save_file($file);
		return $error if $error;
		is_changed("");
	} else {
		my $error = Save_as_file();
		return $error if $error;
	}

	return;
}

sub Save_as_file
{
	if( not is_open() ) {
		return;
	}

	my $file = act_file();
	$file = get_name($file);
	if( $file ) {
		my $error = save_file($file);
		return $error if $error;
	} else {
		my $error = "No file name given.";
		return $error;
	}

	is_saved("YES");
	is_changed("");
	act_file($file);

	return;
}

#------------------------------------------------------------------

sub is_open
{
	return act_file();
}

sub is_changed
{
	return 1;	#dummy
}

sub is_saved
{
	return set_and_get_global("is_saved",$_[0]);

	#my $ref = \$::appl{is_saved};
        #$$ref = $_[0] if defined $_[0];
        #return $$ref;
}

#---------------------

sub act_file
{
	return set_and_get_global("act_grid_file",$_[0]);

	#my $ref = \$::appl{act_grid_file};
        #$$ref = $_[0] if defined $_[0];
        #return $$ref;
}

#---------------------

sub reset_file
{
	Unplot_all();
	act_file("");
}

sub init_file
{
	my $file = shift;

	make_Viewport();
	Plot_all();
	act_file($file) if $file;
}

sub exit_file
{
	Plot_Exit();
}

#---------------------

sub Do_you_want_to_save
{
	my $text = "File has changed. Do you want to save?";
	return Dialog_YesNoCancel($text);
}

sub get_name
{
	my $file = shift;

	return choose_file(".grd",$file);
}

#----------------------------------------------------------

sub new_file
{
	set_grid_and_file(new grd,"Unknown.grd");

	return;
}

sub open_file
{
	my $file = shift;

	print STDERR "open_file...\n";
	my $grid = Read_grid($file);
	return "Error reading file $file" unless $grid;

	set_grid_and_file($grid,$file);

	return;
}

sub close_file
{
	set_grid_and_file(undef,undef);

	return;
}

sub save_file
{
	my $file = shift;

	my $error = Write_grid($file);

	return $error;
}

#------------------------------------------------------------------
1;
#------------------------------------------------------------------

#!/usr/bin/perl -w

use strict;

#-------------------------------------------------------------------------

sub Plot_cross
{
	my ($x,$y,$l,%opt) = @_;
	my $canvas = $::appl{canvas};

	$canvas->createLine($x-$l,$y,$x+$l,$y,%opt);
	$canvas->createLine($x,$y-$l,$x,$y+$l,%opt);
}

sub Plot_circle
{
	my ($x,$y,$r,%opt) = @_;
	my $canvas = $::appl{canvas};

	$canvas->createOval($x-$r,$y-$r,$x+$r,$y+$r,%opt);
}

sub Plot_line_segment
{
	my ($x,$y,$x1,$y1,%opt) = @_;
	my $canvas = $::appl{canvas};

	$canvas->createLine($x,$y,$x1,$y1,%opt);
}

#-------------------------------------------------------------------------

sub Read_grid
{
	my $file = shift;

	my $grid = new grd;

	if( $file ) {
	  $grid->readgrd($file);
	  get_scale_from_grid();
	}

	return $grid;
}

sub Write_grid
{
	my $file = shift;

	my $grid = $::appl{grid};
	return unless $grid;

	$file = "out.grd" unless $file;

	subst_scale_in_grid();
	print STDERR "writing output file...";
	$grid->writegrd($file);

	return;
}

#-------------------------------------------------------------------------

sub Unplot_all
{
	my $canvas = $::appl{canvas};
	$canvas->delete("all");
}

sub Plot_all
{
	Unplot_all();

	Plot_image();
	Plot_grid();

	Plot_tentative();
}

sub Plot_image
{
        my $canpic = $::appl{canpic};
        my $showpic = $::appl{showpic};
	my $canvas = $::appl{canvas};

	my ($cw,$ch) = get_Canvas_Size();

	#$canvas->delete("all");
        $canvas->delete($canpic) if Tk::Exists($canpic);

        $canpic = $canvas->createImage($cw/2,$ch/2,-image=>$showpic);
        $canvas->createRectangle(1,1,$cw-1,$ch-1);

        $::appl{canpic} = $canpic;
}

sub Plot_grid
{
	my $grid = $::appl{grid};
	$grid->make_used();

	print STDERR "plotting grid...\n";

	foreach my $item (values %{$grid->{nodes}}) {
	  my $used = $item->{used};
	  Plot_node($item) unless $used;
	}

	foreach my $item (values %{$grid->{elems}}) {
	  Plot_elem($item);
	}

	foreach my $item (values %{$grid->{lines}}) {
	  Plot_line($item);
	}
}

sub Plot_single_nodes
{
	my $grid = $::appl{grid};
	$grid->make_used();

	print STDERR "plotting single nodes...\n";

	foreach my $n (@_) {
	  my $item = $grid->get_node($n);
	  my $used = $item->{used};
	  Plot_node($item) unless $used;
	}
}

sub Replot_node
{
	my ($node) = @_;

	my $grid = $::appl{grid};
	$grid->make_used();
	my $n = $node->{number};

	Unplot_node($node);
	unless( $node->{used} ) {
	  Plot_node($node);
	  return;
	}

	foreach my $item (values %{$grid->{elems}}) {
	  if( $grid->contains_node($item,$n) ) {
	    Unplot_elem($item);
	    Plot_elem($item);
	  }
	}

	foreach my $item (values %{$grid->{lines}}) {
	  if( $grid->contains_node($item,$n) ) {
	    Unplot_line($item);
	    Plot_line($item);
	  }
	}
}

#-------------------------------------------------------------------------

sub Plot_endpoint
{
	my ($n,%opt) = @_;
	my $radius = 3;
	my $grid = $::appl{grid};

	my $node = $grid->get_node($n);
	my $xr = $node->{x};
	my $yr = $node->{y};
	my ($x,$y) = R2S($xr,$yr);

	Plot_circle($x,$y,$radius,%opt);
}

sub Plot_line
{
	my ($item,%opt) = @_;
	return unless $item;
	my $n = $item->{number};
	unless( $opt{-tags} ) {
	  $opt{-tags} = "Line_$n";
	}
	Plot_elem_line($item,"Line",%opt);
	my $endtag = $::appl{line_end};
	if( $endtag ) {
	  Plot_endpoint($item->{vert}->[0],%opt);
	  Plot_endpoint($item->{vert}->[-1],%opt);
	}
}

sub Plot_elem
{
	my ($item,%opt) = @_;
	return unless $item;
	my $n = $item->{number};
	unless( $opt{-tags} ) {
	  $opt{-tags} = "Elem_$n";
	}
	Plot_elem_line($item,"Elem",%opt);
}

sub Plot_elem_line
{
	my ($item,$what,%opt) = @_;

	my $width = 1;
	my $length = 3;
	my $nodetag = $::appl{plot_node};	#plot node

	return unless $item;
	my $n = $item->{number};
	my $nvert = $item->{nvert};

	my ($xnew,$ynew) = get_xy($item,0);
	($xnew,$ynew) = R2S($xnew,$ynew);
	Plot_cross($xnew,$ynew,$length,-width=>$width,%opt) if $nodetag;

	for(my $i=1;$i<$nvert;$i++) {
	  my ($xold,$yold) = ($xnew,$ynew);
	  ($xnew,$ynew) = get_xy($item,$i);
	  ($xnew,$ynew) = R2S($xnew,$ynew);
	  Plot_line_segment($xold,$yold,$xnew,$ynew,%opt);
	  Plot_cross($xnew,$ynew,$length,-width=>$width,%opt) if $nodetag;
	}

	if( $what eq "Elem" ) {
	  my ($xold,$yold) = R2S(get_xy($item,0));
	  Plot_line_segment($xold,$yold,$xnew,$ynew,%opt);
	}
}

sub get_xy
{
	my ($item,$i) = @_;
	my $grid = $::appl{grid};

	my $vert = $item->{vert};
	my $n = $vert->[$i];

	my $node = $grid->get_node($n);
	return ($node->{x},$node->{y});
}

sub Show_tag_list
{
	my $canvas = $::appl{canvas};
	my @all = $canvas->find("all");

	foreach my $id (@all) {
	  my $type = $canvas->type($id);
	  my @tags = $canvas->gettags($id);
	  my $line = format_tags(@tags);
	  print STDERR "Tag list:  $id  $type - $line\n";
	}
}

sub format_tags
{
	my $line = "";
	foreach my $tag (@_) {
	  $line .= "$tag ";
	}
	return $line;
}

sub is_type
{
	my ($id,$mytype) = @_;
	my $canvas = $::appl{canvas};
	my $type = $canvas->type($id);

	if( $type eq $mytype ) {
	  return 1;
	} else {
	  return 0;
	}
}

sub has_tag
{
	my ($id,$mytag) = @_;
	my $canvas = $::appl{canvas};
	my @tags = $canvas->gettags($id);

	foreach my $tag (@tags) {
	  return 1 if $tag eq $mytag;
	}
	return 0;
}

#-------------------------------------------------------------------

sub Unplot_node_list
{
	my $grid = $::appl{grid};

	foreach my $n (@_) {
	  my $item = $grid->get_node($n);
	  Unplot_node($item);
	}
}

sub Unplot_item
{
	my ($item,$what) = @_;

	return unless $item;

	my $canvas = $::appl{canvas};
	my $n = $item->{number};
	my $tag = "${what}_$n";

	while( my $id = $canvas->find("withtag",$tag) ) {
	  $canvas->delete($id);
	}
}

sub Unplot_node { Unplot_item($_[0],"Node"); }
sub Unplot_elem { Unplot_item($_[0],"Elem"); }
sub Unplot_line { Unplot_item($_[0],"Line"); }

sub Unplot_items
{
	my ($items,$what) = @_;

	foreach my $item (@$items) {
	  print STDERR "unplotting items ($what): $item->{number}\n";
	  Unplot_item($item,$what);
	}
}

sub Unplot_nodes { Unplot_items($_[0],"Node"); }
sub Unplot_elems { Unplot_items($_[0],"Elem"); }
sub Unplot_lines { Unplot_items($_[0],"Line"); }

#-------------------------------------------------------------------

sub Plot_node
{
	my ($item,%opt) = @_;

	return unless $item;

	my $width = 3;
	my $length = 5;

	my $n = $item->{number};
	my $xr = $item->{x};
	my $yr = $item->{y};

	my ($x,$y) = R2S($xr,$yr);

	#print STDERR "plotting node: $xr,$yr  ->  $x,$y\n";

	unless( $opt{-tags} ) {
	  $opt{-tags} = "Node_$n";
	}

	if( $opt{-length} ) {
	  $length = $opt{-length};
	  delete $opt{-length};
	}

	Plot_cross($x,$y,$length,-width=>$width,%opt);
}

#--------------------------------------------------------------

sub Plot_tentative
{
	my $mode = tentative_mode();
	my $item = tentative_item();
	my $extra = tentative_extra();

	Unplot_tentative();

	return unless $item;

	if( $mode eq "TENT_NODE" ) {
	  Plot_tentative_node($item);
	} elsif( $mode =~ "TENT_ELEM" ) {
	  Plot_tentative_elem($item,$extra);
	} elsif( $mode =~ "TENT_LINE" ) {
	  Plot_tentative_line($item,$extra);
	}
}

sub Delete_tentative
{
	reset_tentative_item();
	Unplot_tentative();
}

sub Unplot_tentative
{
        my $canvas = $::appl{canvas};
        $canvas->delete("TentativeItem");
}

sub Plot_tentative_node
{
	my ($item) = @_;
	my $color = "red";

	Plot_node($item,-fill=>$color,-tags=>"TentativeItem");
}

sub Plot_tentative_elem
{
	my ($item,$extra) = @_;
	my $color = "red";

	my $nitem;
	$nitem->{number} = 0;
	$nitem->{x} = $extra->[0];
	$nitem->{y} = $extra->[1];
	Plot_elem($item,-fill=>$color,-tags=>"TentativeItem");
	Plot_node($nitem,-fill=>$color,-tags=>"TentativeItem");
}

sub Plot_tentative_line
{
	my ($item,$extra) = @_;
	my $color = "red";

	my $nitem;
	$nitem->{number} = 0;
	$nitem->{x} = $extra->[0];
	$nitem->{y} = $extra->[1];
	Plot_line($item,-fill=>$color,-tags=>"TentativeItem");
	Plot_node($nitem,-fill=>$color,-tags=>"TentativeItem");
}

sub Delete_act
{
	Delete_tentative();
	Unplot_act();
}

sub Unplot_act
{
        my $canvas = $::appl{canvas};
        $canvas->delete("ActItem");
}

#--------------------------------------------------------------

sub Plot_Postscript
{
        my $canvas = $::appl{canvas};
        print STDERR "plotting Postscript...\n";
        #$canvas->postscript( -file => "out.ps" );
        print STDERR "plotting Postscript finished...\n";
}

sub Plot_Exit
{
	print STDERR "Finished\n";
	my $mw = $::appl{mw};
	$mw->destroy();
}

#--------------------------------------------------------------

sub move_view
{
	my ($xc,$yc,$zoom) = @_;

	if( has_image() ) {
	  move_Image($xc,$yc,$zoom);
	} else {
	  set_new_screen_coords($xc,$yc,$zoom);
	}
}

sub make_Viewport
{
	if( has_image() ) {
	  make_Image();
	} else {
	  my ($xmin,$ymin,$xmax,$ymax) = get_grid_min_max();
	  set_Viewport($xmin,$ymin,$xmax,$ymax);
	  $::appl{original_viewport} = [$xmin,$ymin,$xmax,$ymax];
	}
}

sub set_new_Viewport
{
	if( not has_image() ) {
	  my ($xmin,$ymin,$xmax,$ymax) = get_grid_min_max();
	  set_Viewport($xmin,$ymin,$xmax,$ymax);
	}
}

sub set_original_Viewport
{
	if( not has_image() ) {
	  my ($xmin,$ymin,$xmax,$ymax) = @{$::appl{original_viewport}};
	  set_Viewport($xmin,$ymin,$xmax,$ymax);
	}
}

sub set_Viewport
{
	my ($xmin,$ymin,$xmax,$ymax) = @_;

	if( not has_image() ) {
	  my $scale = 0;
	  ($xmin,$ymin,$xmax,$ymax,$scale) = 
		make_Viewport_extremes($xmin,$ymin,$xmax,$ymax);
	  set_scale(0,0,$xmin,$ymax,$scale);
	}
}

sub make_Viewport_extremes
{
	my ($xmin,$ymin,$xmax,$ymax) = @_;

	my $fac = 0.1;

	my $dx = $xmax-$xmin;
	my $dy = $ymax-$ymin;
	my $dxy = max($dx,$dy);
	my $scale;
	my ($px,$py) = get_Viewport_Size();

	if( $dx/$px > $dy/$py ) {
	  $scale = $dx*(1+2*$fac)/$px;
	  $ymin -= 0.5 * ($py*$dx/$px-$dy);
	} else {
	  $scale = $dy*(1+2*$fac)/$py;
	  $xmin -= 0.5 * ($px*$dy/$py-$dx);
	}

	$xmin -= $fac*$dxy;
	$xmax += $fac*$dxy;
	$ymin -= $fac*$dxy;
	$ymax += $fac*$dxy;

	return ($xmin,$ymin,$xmax,$ymax,$scale);
}

#--------------------------------------------------------------

sub write_coords
{
	my ($xr,$yr) = @_;

	my ($xs,$ys) = R2S($xr,$yr);
	my ($xp,$yp) = R2P($xr,$yr);

	print "coords (r,p,s): $xr,$yr   $xp,$yp   $xs,$ys\n";
}

sub write_hash
{
	my %hash = @_;

	foreach my $key (keys %hash) {
	  print STDERR "$key : $hash{$key}\n";
	}
}

#----------------------
1;
#----------------------

#!/usr/bin/perl -w

use strict;

#--------------------------------------------------------------------

sub make_node
{
	my ($type,$x,$y) = @_;

	my $item;

	my $grid = $::appl{grid};
	my $depth_on = get_global("depth_on");
	my $depthvalue = get_global("depthvalue");

	if( $depth_on ) {
	  $item = $grid->make_node(0,$type,$x,$y,$depthvalue);
	} else {
	  $item = $grid->make_node(0,$type,$x,$y);
	}

	return $item;
}

sub delete_node
{
	my ($item) = @_;

	my $grid = $::appl{grid};

	$grid->delete_node($item);
}

sub unify_nodes
{
	my ($item1,$item2) = @_;

	my $grid = $::appl{grid};

	$grid->unify_nodes($item1,$item2);
}

sub is_used
{
	my ($node) = @_;

	my $grid = $::appl{grid};
	$grid->make_used();

        return $node->{used};
}

#--------------------------------------------------------------------

sub make_elem
{
	my ($rf) = @_;

	my $item;

	my @f = @$rf;
	my $n = @f;
	return if $n < 3;

	my @numbs = ();
	foreach $item (@f) {
	  push(@numbs,$item->{number});
	}

	my $grid = $::appl{grid};
	my $depth_on = get_global("depth_on");
	my $depthvalue = get_global("depthvalue");

	if( $depth_on ) {
	  $item = $grid->make_elem(0,0,$depthvalue,@numbs);
	} else {
	  $item = $grid->make_elem(0,0,"",@numbs);
	}
	print STDERR "writing elem: $n  ($item->{number})\n";
	my $line = join(" ",@numbs);
	print STDERR "$line\n";

	return $item;
}

#--------------------------------------------------------------------

sub make_line
{
	my ($rf) = @_;

	my $item;

	my @f = @$rf;
	my $n = @f;
	return if $n < 2;
	pop(@f) if $f[$n-1] eq $f[$n-2];
	$n = @f;
	return if $n < 2;

	my @numbs = ();
	foreach $item (@f) {
	  push(@numbs,$item->{number});
	}

	my $grid = $::appl{grid};
	my $depth_on = get_global("depth_on");
	my $depthvalue = get_global("depthvalue");

	if( $depth_on ) {
	  $item = $grid->make_line(0,0,$depthvalue,@numbs);
	} else {
	  $item = $grid->make_line(0,0,"",@numbs);
	}
	print STDERR "writing line: $n  ($item->{number})\n";
	my $line = join(" ",@numbs);
	print STDERR "$line\n";

	return $item;
}

sub delete_lines
{
	my ($lines) = @_;

	my $grid = $::appl{grid};

	foreach my $line (@$lines) {
	  $grid->delete_line($line);
	}
}

sub delete_line
{
	my ($item) = @_;

	my $grid = $::appl{grid};

	$grid->delete_line($item);
}

sub close_line
{
	my ($item) = @_;

	return unless $item;
	my $grid = $::appl{grid};

	$grid->close_line($item);
}

sub auto_close_line
{
	my ($item) = @_;

	close_line($item) if get_global("auto_close");
}

sub split_line
{
	my ($item,$node) = @_;

	my $grid = $::appl{grid};

	return $grid->split_line($item,$node);
}

sub join_line
{
	my ($line1,$line2,$node) = @_;

	my $grid = $::appl{grid};

	return $grid->connect_lines($line1,$line2,$node);
}

sub join_all
{
	my ($line) = @_;

	my $v = $line->{vert};
	my $l1 = [];
	my $n1 = [];
	my ($l2,$n2,$err);
	my $closed = 0;

	print STDERR "join all 1\n";
	($l2,$n2,$err) = get_lines_from_node($v->[-1],$v->[0],$line);
	return [] if $err;
	$closed = 1 if scalar(@$l2) and $n2->[-1] == $v->[0];
	print STDERR "join all 2 ($err $closed)\n";
	unless( $closed ) {
	  ($l1,$n1,$err) = get_lines_from_node($v->[0],$v->[-1],$line);
	  print STDERR "join all 3 ($err)\n";
	  return [] if $err;
	}

	print_array("first line nodes:  ",@$n2);
	print_array("second line nodes: ",@$n1);

	my @lines = ();
	push(@lines,@$l1);
	push(@lines,@$l2);

	my @nodes = ();
	push(@nodes,reverse(@$n1));
	push(@nodes,@$v);
	push(@nodes,@$n2);

	print_array("all nodes:  ",@nodes);
	@nodes = elim_double_nodes(@nodes);
	print_array("elim nodes:  ",@nodes);
	my $n = @lines;
	print STDERR "$n lines to be deleted\n";

	$line->{vert} = \@nodes;
	$line->{nvert} = scalar @nodes;

	return \@lines;		#still to delete
}

sub get_lines_from_node
{
	my ($start,$stop,$line) = @_;

	my $grid = $::appl{grid};
	my $nline = $line->{number};
	my $laux = $line;		#last line inserted
	my @nodes = ();
	my @lines = ();
	my $count;

	do {

	my $nlast = $laux->{number};
	$count = 0;
	my @naux = ();

	print STDERR "join: $nlast $nline\n";

        foreach my $item (values %{$grid->{lines}}) {
          my $n = $item->{number};
	  next if $n == $nline or $n == $nlast;
          my $v = $item->{vert};
	  if( $v->[0] == $start or $v->[-1] == $start ) {
	    @naux = @$v;
	    @naux = reverse @$v if $v->[-1] == $start;
	    $laux = $item;
	    $count++;
	  }
	}

	if( $count > 1 ) {
	  return ([],[],-1);	#error
	} elsif( $count == 1 ) {
	  push(@nodes,@naux);
	  push(@lines,$laux);
	  $start = $naux[-1];
	}

	print STDERR "   join: $count $laux->{number}\n";

	} while( $count and $nodes[-1] != $stop );

	return (\@lines,\@nodes,0);
}

sub elim_double_nodes
{
	my @new = ($_[0]);

	foreach my $n (@_) {
	  push(@new,$n) unless $new[-1] == $n;
	}

	return @new;
}

sub delete_single_nodes
{
        my $grid = $::appl{grid};
	$grid->make_used();

	foreach my $n (@_) {
          my $item = $grid->get_node($n);
          my $used = $item->{used};
	  $grid->delete_node($item) unless $used;
	}
}

#-------------------------------------------------------------

sub get_closest_node
{
        my ($xr,$yr,$old) = @_;

        my $grid = $::appl{grid};

        my $dmin = 1.e30;
	my @item_list = ();

        foreach my $item (values %{$grid->{nodes}}) {
          my ($x,$y) = ($item->{x},$item->{y});
          my $dist = ($x-$xr)*($x-$xr) + ($y-$yr)*($y-$yr);
          if( $dist < $dmin ) {
            $dmin = $dist;
	    @item_list = ($item);
          } elsif( $dist == $dmin ) {
	    push(@item_list,$item);
          }
        }

	return get_next_item($old,@item_list);
}

sub get_next_item
{
	my $old = shift;
	return $_[0] unless $old;

	my $ntot = @_;

	if( $ntot == 1 ) {		#only one item
	  return $_[0];
	} elsif( $ntot == 0 ) {		#no item
	  return;
	}
	  
	for(my $i=0;$i<$ntot;$i++) {
	  if( $_[$i] eq $old ) {
	    my $j = ($i+1) % $ntot;
	    return $_[$j];		#next item
	  }
	}

	return $_[0];			#not equal to old item
}

#-------------------------------------------------------

sub get_closest_line
{
        my ($xr,$yr,$old) = @_;

        my $grid = $::appl{grid};
	my $mode = tentative_mode();
	$mode =~ s/TENT_LINE//;

	my @item_list = ();
        my $dmin = 1.e30;
	#my $d;
	my ($xc,$yc,$nmin,$smin);
	my ($xn,$yn,$node,$seg);

        foreach my $item (values %{$grid->{lines}}) {
          my $v = $item->{vert};
	  my ($d,$xn,$yn) = get_closest_point_on_line($xr,$yr,$v);
          if( $d < $dmin ) {
            $dmin = $d;
	    ($xc,$yc) = ($xn,$yn);
	    @item_list = ($item);
          } elsif( $d == $dmin ) {
	    push(@item_list,$item);
          }
	}

	my $imin = get_next_item($old,@item_list);
	return unless $imin;
        my $v = $imin->{vert};
	if( $mode eq "_NODE" ) {
	    ($xc,$yc,$node) = find_closest_node($xc,$yc,$v);
	} elsif( $mode eq "_SEGMENT" ) {
	    ($xc,$yc,$seg) = find_closest_segment($xc,$yc,$v);
	}

	return ($imin,$xc,$yc,$node,$seg);
}

sub find_closest_node
{
	my ($x,$y,$rl) = @_;
	my $grid = $::appl{grid};

        my $dmin = 1.e30;
	my ($xc,$yc);
	my $imin = "";
	my $n = @$rl;

	for( my $i=0;$i<$n;$i++ ) {
	  my $node = $grid->get_node($rl->[$i]);  
	  my ($xn,$yn) = ($node->{x},$node->{y});
	  my $d = ($xn-$x)*($xn-$x) + ($yn-$y)*($yn-$y);
	  if( $d < $dmin ) {
	    $dmin = $d;
	    ($xc,$yc) = ($xn,$yn);
	    $imin = $node;
	  }
	}
	return ($xc,$yc,$imin);
}

sub find_closest_segment
{
	my ($x,$y,$rl) = @_;
	my $grid = $::appl{grid};

        my $dmin = 1.e30;
	my ($xc,$yc);
	my $iseg = 0;
	my $n = @$rl;

	for( my $i=1;$i<$n;$i++ ) {
	  my $node1 = $grid->get_node($rl->[$i-1]);  
	  my ($x1,$y1) = ($node1->{x},$node1->{y});
	  my $node2 = $grid->get_node($rl->[$i]);  
	  my ($x2,$y2) = ($node2->{x},$node2->{y});
	  my ($xn,$yn) = (0.5*($x1+$x2),0.5*($y1+$y2));
	  my $d = ($xn-$x)*($xn-$x) + ($yn-$y)*($yn-$y);
	  if( $d < $dmin ) {
	    $dmin = $d;
	    ($xc,$yc) = ($xn,$yn);
	    $iseg = $i;
	  }
	}
	return ($xc,$yc,$iseg);
}

sub get_closest_point_on_line
{
	my ($x,$y,$rl) = @_;
	my $grid = $::appl{grid};

        my $dmin = 1.e30;
	my ($xc,$yc);
	my $n = @$rl;

	for( my $i=1;$i<$n;$i++ ) {
	  my $node1 = $grid->get_node($rl->[$i-1]);  
	  my ($x1,$y1) = ($node1->{x},$node1->{y});
	  my $node2 = $grid->get_node($rl->[$i]);  
	  my ($x2,$y2) = ($node2->{x},$node2->{y});
	  my ($xn,$yn,$s) = get_closest_point_on_segment($x,$y,$x1,$y1,$x2,$y2);
	  my $d = ($xn-$x)*($xn-$x) + ($yn-$y)*($yn-$y);
	  if( $d < $dmin ) {
	    $dmin = $d;
	    ($xc,$yc) = ($xn,$yn);
	  }
	}

	return ($dmin,$xc,$yc);
}

sub get_closest_point_on_segment
{
	my ($x,$y,$x1,$y1,$x2,$y2) = @_;

	my $x12 = $x2 - $x1;
	my $y12 = $y2 - $y1;
	my $x1p = $x - $x1;
	my $y1p = $y - $y1;

	my $d12 = $x12*$x12+$y12*$y12;
	my $s = 0.;
	if( $d12 > 0 ) {
	  $s = ($x1p*$x12+$y1p*$y12)/$d12;
	}

	$s = 0 if $s < 0;
	$s = 1 if $s > 1;

	$x = $x1 + $s * $x12;
	$y = $y1 + $s * $y12;

	return ($x,$y,$s);
}

#-------------------------------------------------------------

sub get_closest_elem
{
        my ($xr,$yr,$old) = @_;

        my $grid = $::appl{grid};
	my $mode = tentative_mode();
	$mode =~ s/TENT_ELEM//;

	my @item_list = ();
        my $dmin = 1.e30;
	#my $d;
	my ($xc,$yc,$nmin,$smin);
	my ($xn,$yn,$node,$seg);

        foreach my $item (values %{$grid->{elems}}) {
	  if( point_in_element($xr,$yr,$item) ) {
	    push(@item_list,$item);
          }
	}

	my $imin = get_next_item($old,@item_list);

	$grid->elem_info($old,"old elem: ");
	$grid->elem_info($imin,"new elem: ");
	return ("") unless $imin;

        my $v = $imin->{vert};
	if( $mode eq "_NODE" ) {
	    ($xc,$yc,$node) = find_closest_node($xc,$yc,$v);
	} elsif( $mode eq "_SEGMENT" ) {
	    my @vaux = @$v;
	    push(@vaux,$vaux[0]);	#attach first point
	    ($xc,$yc,$seg) = find_closest_segment($xc,$yc,\@vaux);
	} else {
	    ($xc,$yc) = make_center_point($v);
	}

	return ($imin,$xc,$yc,$node,$seg);
}

sub make_center_point
{
	my $nlist = shift;

        my $grid = $::appl{grid};
	my ($xc,$yc,$nn) = (0,0,0);

	foreach my $n (@$nlist) {
	  my $item = $grid->get_node($n);
	  $xc += $item->{x};
	  $yc += $item->{y};
	  $nn++;
	}

	if( $nn ) {
	  $xc /= $nn;
	  $yc /= $nn;
	}

	return ($xc,$yc);
}

sub point_in_element
{
	my ($xr,$yr,$item) = @_;

        my $grid = $::appl{grid};
	my ($x,$y) = $grid->make_xy($item);
	my $gl = new grdline;
	$gl->set_line($x,$y);

	my $return =  $gl->in_line($xr,$yr);

	return $return;
}

sub delete_elem
{
	my ($item) = @_;

	my $grid = $::appl{grid};

	$grid->delete_elem($item);
}

#-------------------------------------------------------------

sub get_grid_min_max
{
	my $grid = $::appl{grid};
	my @items = values %{$grid->{nodes}};

	my ($xmin,$ymin,$xmax,$ymax);

	my $item = $items[0];
	if( $item ) {
	  ($xmin,$ymin) = ($item->{x},$item->{y});
	  ($xmax,$ymax) = ($xmin,$ymin);
	} else {
	  ($xmin,$ymin,$xmax,$ymax) = (0,0,1,1);
	}

        foreach my $item (@items) {
	  my ($x,$y) = ($item->{x},$item->{y});
	  $xmin = $x if $x < $xmin;
	  $ymin = $y if $y < $ymin;
	  $xmax = $x if $x > $xmax;
	  $ymax = $y if $y > $ymax;
	}

	if( $xmin == $xmax ) { ($xmin,$xmax) = ($xmin-1,$xmax+1); }
	if( $ymin == $ymax ) { ($ymin,$ymax) = ($ymin-1,$ymax+1); }

	return ($xmin,$ymin,$xmax,$ymax);
}

#-------------------------------------------------------------

sub rescale_points
{
	my $scale = shift;

        my $grid = $::appl{grid};
        my $nodes = $grid->{nodes};

        foreach my $item (values %$nodes) {
          my $x = $item->{x};
          my $y = $item->{y};
          my ($xn,$yn) = &$scale($x,$y);
          $item->{x} = $xn;
          $item->{y} = $yn;
          print STDERR "scaling  $x $y  ->  $xn $yn\n";
        }
}

#-------------------------------------------------------------

sub subst_scale_in_grid
{
	my $grid = $::appl{grid};
	my $comments = $grid->{comms};

	my @new = ();

	foreach my $line (@$comments) {		#eliminate all SCALE directives
	  unless( $line =~ /^0 \(SCALE\)\s+/ ) {
	    push(@new,$line);
	  }
	}

	my ($xp0,$yp0,$xr0,$yr0,$scale) = get_scale();

	my $line = "0 (SCALE) $xp0 $yp0 $xr0 $yr0 $scale\n";
	push(@new,$line);

	$grid->{comms} = \@new;
}

sub get_scale_from_grid
{
	my $grid = $::appl{grid};
	my $comments = $grid->{comms};

	foreach my $line (@$comments) {
	  if( $line =~ /^0 \(SCALE\)\s+(.*)$/ ) {
	    my @f = split(/\s+/,$1);
	    if( scalar(@f) != 5 ) {
		die "cannot read scale: $line\n";
	    }
	    set_scale(@f);
	  }
	}
}

#-------------------------------------------------------------

#----------------------
1;
#----------------------

#!/usr/bin/perl -w

use strict;

#---------------------------------------------------------------

sub toggle
{
	return ($_[0]+1)%2;
}

sub max
{
	my $max = shift;

	foreach my $i (@_) {
	  $max = $i if $i > $max;
	}

	return $max;
}

sub min
{
	my $min = shift;

	foreach my $i (@_) {
	  $min = $i if $i < $min;
	}

	return $min;
}

sub dist
{
	my ($x1,$y1,$x2,$y2) = @_;

	my $dx = $x1-$x2;
	my $dy = $y1-$y2;

	my $dist = $dx*$dx + $dy*$dy;

	return sqrt($dist);
}

sub round_binary
{
	my ($frac) = @_;

	my $zoom = 1;

	if( $frac < 1 ) {
	  while( $zoom > $frac ) {
	    $zoom /= 2;
	  }
	} elsif( $frac > 1 ) {
	  while( 2*$zoom < $frac ) {
	    $zoom *= 2;
	  }
	} else {
	  ;
	}

	return $zoom;
}

sub print_array
{
	my $text = shift;

	print "$text\n";
	foreach my $item (@_) {
	  print "$item ";
	}
	print "\n";
}

#--------------------------------------------------------------

sub set_scale
{
	my ($xp0,$yp0,$xr0,$yr0,$scale) = @_;

	$::appl{scale_info} = [$xp0,$yp0,$xr0,$yr0,$scale];
}

sub get_scale
{
	return @{$::appl{scale_info}};
}

#--------------------------------------------------------------

sub R2P		# transforms real coordinates to picture coordinates
{
	my ($xr,$yr) = @_;

	my ($xp0,$yp0,$xr0,$yr0,$scale) = @{$::appl{scale_info}};

	my $xp = ($xr-$xr0)/$scale + $xp0;
	my $yp = -($yr-$yr0)/$scale + $yp0;

	return ($xp,$yp);
}

sub P2R		# transforms picture coordinates to screen coordinates
{
	my ($xp,$yp) = @_;

	my ($xp0,$yp0,$xr0,$yr0,$scale) = @{$::appl{scale_info}};

	my $xr = ($xp-$xp0)*$scale + $xr0;
	my $yr = -($yp-$yp0)*$scale + $yr0;

	return ($xr,$yr);
}

sub R2S
{
	my ($xr,$yr) = @_;
	return P2S(R2P($xr,$yr));
}

sub S2R
{
	my ($xs,$ys) = @_;
	return P2R(S2P($xs,$ys));
}

sub Rescale_points
{

	my ($inverse) = @_;
	my $scale;

	if( $inverse ) {
	  $scale = \&R2P;
	} else {
	  $scale = \&P2R;
	}

	rescale_points($scale);
}

#-------------------------------------------------------------

#----------------------
1;
#----------------------

#!/usr/bin/perl -w
#
# grd utility routines
#
##############################################################
#
# version 1.4.1
#
# 19.08.2005		unify_nodes, if defined $depth
# 24.08.2005		connect_lines, split_line, contains_node
# 26.08.2005		make_xy, node_info, elem_info
# 20.10.2005		version control introduced
#
##############################################################

# example of usage:
#
# use lib "$ENV{HOME}/shyfem/femlib/perl";
# use grd;
# 
# my $grid = new grd;
# my $infile = $ARGV[0];
# 
# $grid->readgrd($infile);
# $grid->writegrd("out.grd");
# 
# to iterate over say the list of elements:
#
# my $elems = $grid->get_elems();
# foreach my $item (values %$elems) {
#   my $number = $item->{number};
#   ...
# }
#
# or
#
# foreach my $item ( $grid->get_elem_list() ) {
#   my $number = $item->{number};
#   ...
# }
#
##############################################################

use strict;

package grd;

##############################################################

sub new
{
    my $self;

    my %nodes = ();
    my %elems = ();
    my %lines = ();
    my @comms = ();

    $self =	{
	    		 file		=>	undef
	    		,outfile	=>	undef
	    		,outhandle	=>	*STDOUT
			,nodes		=>	\%nodes
			,elems		=>	\%elems
			,lines		=>	\%lines
			,comms		=>	\@comms
			,nnmax		=>	0
			,nemax		=>	0
			,nlmax		=>	0
			,verbose	=>	1
		};

    bless $self;
    return $self;
}

###############################################################################

sub set_outfile {

  my ($self,$outfile) = @_;

  if( $self->{outhandle} and $self->{outhandle} ne *STDOUT ) {
      close($self->{outhandle});
  }

  my $handle;
  if( $outfile ) {
    my $file = make_grdfile($outfile);
    open($handle,">$outfile") || die "Cannot open outfile: $file\n";
  } else {
    $handle = *STDOUT;
  }

  $self->{outfile} = $outfile;
  $self->{outhandle} = $handle;
}

sub writegrd {

  my ($self,$file) = @_;

  $self->set_outfile($file);
  my $fh = $self->{outhandle};
  $file = $self->{outfile};

  if( $self->{verbose} ) {
    print STDERR "writing file: $file\n" if $file;
  }

  my ($items,@keys);

  print $fh "\n";

  foreach my $line (@{$self->{comms}}) {
    print $fh "$line\n";
  }
  print $fh "\n";
    
  $items = $self->{nodes};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  foreach my $number (@keys) {
      $self->write_node($items->{$number});
  } 
  print $fh "\n";

  $items = $self->{elems};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  foreach my $number (@keys) {
      $self->write_elem($items->{$number});
  } 
  print $fh "\n";

  $items = $self->{lines};
  @keys = keys %$items;
  @keys = sort {$a<=>$b} @keys;
  foreach my $number (@keys) {
      $self->write_line($items->{$number});
  } 
  print $fh "\n";

}

sub write_nodes {		# helper routine

  my ($fh,$ra) = @_;
  my $j = 0;

  foreach (@$ra) {
    $j++;
    print $fh " $_";
    print $fh "\n" if $j%10 == 0;
  }
}

sub print_node
{
    my ($self,$number) = @_;

    my $node = $self->{nodes}->{$number};

    $self->write_node($node);
}

sub write_node
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};

    print $fh "1 $item->{number} $item->{type}";
    print $fh " $item->{x} $item->{y}";
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n";
}

sub write_elem
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};
    my $nvert = $item->{nvert};

    print $fh "2 $item->{number} $item->{type}";
    print $fh " $item->{nvert}";
    print $fh "\n" if $nvert > 3;
    &write_nodes($fh,$item->{vert});
    print $fh "\n" if $nvert > 3 and $item->{h};
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n";
} 

sub write_line
{
    my ($self,$item) = @_;

    my $fh = $self->{outhandle};

    print $fh "3 $item->{number} $item->{type} $item->{nvert}\n";
    &write_nodes($fh,$item->{vert});
    print $fh " $item->{h}" if defined $item->{h};
    print $fh "\n\n";
}

###############################################################################

sub readgrd {

  my ($self,$file) = @_;

  return unless $file;

  $file = make_grdfile($file);
  open(FILE,"$file") || die "Cannot open file: $file\n";

  $self->{file} = $file;
  if( $self->{verbose} ) {
    print STDERR "reading file: $file\n";;
  }

  while( $_ = $self->nextitem ) {

	#print "readgrd: $_\n";

        next if /^\s*$/;
	my @f = split;

	my $item = $f[0];

	if( $item == 0 ) {
		$self->insert_comment($_);
	} elsif( $item == 1 ) {
		$self->insert_node(\@f);
	} elsif( $item == 2 ) {
		$self->insert_elem(\@f);
	} elsif( $item == 3 ) {
		$self->insert_line(\@f);
	} else {
		die "Unknown item: $_\n";
	}
  }

  close(FILE);
}

sub make_grdfile {

  my $name = shift;
  my $file = "";

  if( $name =~ /\.grd\s*$/i ) {
    $file = $name;
  } else {
    $file = "$name.grd";
  }

  return $file;
}

###################################

sub nextitem {

  my ($self) = @_;

  my $line = &getline;
  return undef unless defined($line);

  #print "nextitem (line): $line";

  while( my $newline = &getline ) {
    #print "nextitem (newline): $newline";
    if( $newline =~ /^\s*$/ ) {	#empty line
	last;
    } elsif( $newline =~ /^\s+/ ) {		#conti line
	$line .= " $newline";
    } else {
	&ungetline($newline);		#save for next call
	last;
    }
  }

  $line =~ s/\n/ /g;

  return $line;
}

###################################

sub getline {

  if( @grd::oldlines ) {
    return shift(@grd::oldlines);
  } else {
    return <FILE>;
  }
}

sub ungetline {

  push(@grd::oldlines,$_[0]);
}

###################################

sub insert_comment 
{
    my ($self,$line) = @_;

    my $rcomms = $self->{comms};
    push(@$rcomms,$line);
}

sub insert_node 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $rnodes = $self->{nodes};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{x}      = $$ritems[3];
    $item{y}      = $$ritems[4];
    $item{h}      = $$ritems[5];

    $item{used}   = 0;
    $self->{nnmax} = $number if $number > $self->{nnmax};

    $rnodes->{$number} = \%item;
}

sub insert_elem 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $nvert  = $$ritems[3];
    my $relems = $self->{elems};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{nvert}  = $nvert;
    $item{vert}   = &get_vertices($ritems);
    $item{h}      = $$ritems[4+$nvert];

    $self->{nemax} = $number if $number > $self->{nemax};

    $relems->{$number} = \%item;
}

sub insert_line 
{
    my ($self,$ritems) = @_;

    my %item = ();

    my $number = $$ritems[1];
    my $nvert  = $$ritems[3];
    my $rlines = $self->{lines};

    $item{number} = $number;
    $item{type}   = $$ritems[2];
    $item{nvert}  = $nvert;
    $item{vert}   = &get_vertices($ritems);
    $item{h}      = $$ritems[4+$nvert];

    $self->{nlmax} = $number if $number > $self->{nlmax};

    $rlines->{$number} = \%item;
}

###############################################################################

sub node_info
{
    my ($self,$item,$text) = @_;

    print STDERR "$text\n";
    if( $item ) {
      print STDERR "$item->{number} $item->{type} $item->{x} $item->{y}\n";
    } else {
      print STDERR "no such item\n";
    }
}

sub elem_info
{
    my ($self,$item,$text) = @_;

    print STDERR "$text\n";
    if( $item ) {
      print STDERR "$item->{number} $item->{type} $item->{nvert}\n";
    } else {
      print STDERR "no such item\n";
    }
}

###############################################################################

sub by_number {
  $a<=>$b;
}

sub compress_node_numbers
{
    my ($self) = @_;

    my $nodes = $self->{nodes};
    my @numbers = sort by_number keys(%$nodes);
    my %intern = ();
    my %new_nodes = ();
    my $items;

    my $i = 1;
    foreach my $number (@numbers) {
      $intern{$number} = $i;
      my $node = $nodes->{$number};
      $node->{number} = $i;
      $new_nodes{$i} = $node;
      $i++;
    }
    $self->{nodes} = \%new_nodes;

    $items = $self->{elems};
    &compress_node_list($items,\%intern);

    $items = $self->{lines};
    &compress_node_list($items,\%intern);
}

sub compress_node_list
{

    my ($items,$intern) = @_;

    my @keys = keys %$items;
    foreach my $number (@keys) {
      my $item = $items->{$number};
      my $node_list = $item->{vert};
      my @new_nodes = ();
      foreach my $node (@$node_list) {
	push(@new_nodes,$$intern{$node});
      }
      $item->{vert} = \@new_nodes;
    }
}

###############################################################################

sub get_node_list { return $_[0]->get_item_list("nodes"); }
sub get_elem_list { return $_[0]->get_item_list("elems"); }
sub get_line_list { return $_[0]->get_item_list("lines"); }

sub get_item_list
{
    my ($self,$type) = @_;

    my $items = $self->{$type};
    return values %$items
}

#----------

sub get_nodes { return $_[0]->get_items("nodes"); }
sub get_elems { return $_[0]->get_items("elems"); }
sub get_lines { return $_[0]->get_items("lines"); }

sub get_items
{
    my ($self,$type) = @_;

    return $self->{$type};
}

#----------

sub exists_node { return $_[0]->exists_item($_[1],"nodes"); }
sub exists_elem { return $_[0]->exists_item($_[1],"elems"); }
sub exists_line { return $_[0]->exists_item($_[1],"lines"); }

sub exists_item
{
    my ($self,$number,$type) = @_;

    my $items = $self->{$type};
    if( exists $items->{$number} ) {
      return 1;
    } else {
      return 0;
    }
}

#----------

sub get_node { return $_[0]->get_item($_[1],"nodes"); }
sub get_elem { return $_[0]->get_item($_[1],"elems"); }
sub get_line { return $_[0]->get_item($_[1],"lines"); }

sub get_item
{
    my ($self,$number,$type) = @_;

    my $items = $self->{$type};
    unless( exists $items->{$number} ) {
      die "get_item: Cannot find number $number of type $type\n";
    }
    return $items->{$number};
}

#----------
#
# $grid->delete_node($node)    where $node is ref to node item
# $grid->delete_node(511)      directly with number

sub delete_node { return $_[0]->delete_item($_[1],"nodes"); }
sub delete_elem { return $_[0]->delete_item($_[1],"elems"); }
sub delete_line { return $_[0]->delete_item($_[1],"lines"); }

sub delete_item
{
    my ($self,$item,$type) = @_;

    my $items = $self->{$type};
    $item = $self->get_item($item,$type) unless ref($item); #if number, get item
    my $number = $item->{number};

    delete $$items{$number};
}

#----------

sub clone_node { return $_[0]->clone_item($_[1],"nodes"); }
sub clone_elem { return $_[0]->clone_item($_[1],"elems"); }
sub clone_line { return $_[0]->clone_item($_[1],"lines"); }

sub clone_item
{
    my ($self,$item,$type) = @_;

    my $items = $self->{$type};
    my $number = $item->{number};
    if( exists $items->{$number} ) {
      print STDERR "clone_item: Item $number of type $type exists already\n";
    } else {
      $items->{$number} = $item;
    }
}

###############################################################################

sub unify_nodes {	#unifies nodes -> node n2 is deleted

  my ($self,$n1,$n2) = @_;

  $n1 = $n1->{number} if ref($n1);	# get number if item
  $n2 = $n2->{number} if ref($n2);	# get number if item

  $self->delete_node($n2);
  $self->substitute_node($n1,$n2);
}

sub substitute_node {

  my ($self,$n1,$n2) = @_;

  # substitute n2 with n1

  $self->substitute_vert($self->{elems},$n1,$n2);
  $self->substitute_vert($self->{lines},$n1,$n2);
}

sub substitute_vert {

  my ($self,$items,$n1,$n2) = @_;

  # substitute n2 with n1

  foreach my $item (values %$items) {
    my $vert = $item->{vert};
    foreach my $n (@$vert) {
	$n = $n1 if $n == $n2;
    }
  }
}

###################################

sub make_node
{
    my ($self,$number,$type,$x,$y,$depth) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nnmax}++;
	$number = $self->{nnmax};
    }

    $item{number} = $number;
    $item{type}   = $type;
    $item{x}      = $x;
    $item{y}      = $y;
    $item{h}      = $depth if defined $depth;

    $self->{nodes}->{$number} = \%item;

    return \%item;
}

sub make_elem
{
    my ($self,$number,$type,$depth,@vert) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nemax}++;
	$number = $self->{nemax};
    }

    #print STDERR "make_elem: $vert[0] $vert[1] $vert[2]\n";

    $item{number} = $number;
    $item{type}   = $type;
    $item{nvert}  = @vert;
    $item{vert}   = \@vert;
    $item{h}      = $depth if defined $depth;

    $self->{elems}->{$number} = \%item;

    return \%item;
}

sub make_line
{
    my ($self,$number,$type,$depth,@vert) = @_;

    my %item = ();

    if( $number <= 0 ) {
	$self->{nlmax}++;
	$number = $self->{nlmax};
    }

    $item{number} = $number;
    $item{type}   = $type;
    $item{nvert}  = @vert;
    $item{vert}   = \@vert;
    $item{h}      = $depth if defined $depth;

    $self->{lines}->{$number} = \%item;

    return \%item;
}

###################################

sub make_xy
{
	my ($self,$item) = @_;

	my @x = ();
	my @y = ();
	my $v = $item->{vert};

	foreach my $n (@$v) {
	  my $node = $self->get_node($n);
	  push(@x,$node->{x});
	  push(@y,$node->{y});
	}

	return (\@x,\@y);
}

sub contains_node
{
  my ($self,$item,$n) = @_;

  my $v = $item->{vert};

  foreach my $node (@$v) {
    return 1 if $node == $n;
  }
  return 0;
}

sub is_closed
{
  my ($self,$line) = @_;

  my $vert = $line->{vert};

  if( $$vert[0] == $$vert[-1] ) {
	return 1;
  } else {
	return 0;
  }
}

sub close_line
{
  my ($self,$line) = @_;

  my $vert = $line->{vert};

  if( $$vert[0] != $$vert[-1] ) {
    push(@$vert,$$vert[0]);
    $line->{nvert}++;
  }
}

sub connect_lines
{
  # we try to use node, but if node is wrong we connect anyway

  my ($self,$line1,$line2,$node) = @_;

  my $vert1 = $line1->{vert};
  my $vert2 = $line2->{vert};
  my $n = $node->{number} if $node;

  my @new = ();
  my @v1 = ();
  my @v2 = ();

  if( $$vert1[-1] == $$vert2[0] ) {
    @v1 = @$vert1;
    @v2 = @$vert2;
  } elsif( $$vert1[0] == $$vert2[0] ) {
    @v1 = reverse @$vert1;
    @v2 = @$vert2;
  } elsif( $$vert1[-1] == $$vert2[-1] ) {
    @v1 = @$vert1;
    @v2 = reverse @$vert2;
  } elsif( $$vert1[0] == $$vert2[-1] ) {
    @v1 = reverse @$vert1;
    @v2 = reverse @$vert2;
  } else {
    print STDERR  "*** Cannot connect lines $line1 and $line2: no node in common\n";
    return;
  }

  if( $node and $v2[0] != $n ) {	#node is given and link not ok
    if( $v1[0] == $n and $v2[-1] == $n ) {
	my @aux = @v1;
	@v1 = @v2;
	@v2 = @aux;
    }
  }
  push(@new,@v1);
  shift(@v2);
  push(@new,@v2);

  $line1->{vert} = \@new;
  $line1->{nvert} = @new;

  $self->delete_line($line2);

  return $line1;
}

sub split_line
{
  my ($self,$line,$node) = @_;

  my $number = $node->{number};
  my $vert = $line->{vert};
  my @v1 = ();
  my @v2 = ();
  my $act = \@v1;

  foreach my $n (@$vert) {
	push(@$act,$n);
	if( $n == $number ) {
	  $act = \@v2;
	  push(@$act,$n);
	}
  }
	  
  if( scalar(@v1) < 2 or scalar(@v2) < 2 ) {
    print STDERR "*** cannot split line on this point\n";
    return ($line,"");
  }

  $line->{nvert} = scalar(@v1);
  $line->{vert} = \@v1;

  my $line_new = $self->make_line(0,$line->{type},$line->{depth},@v2);

  return ($line,$line_new);
}

###################################

sub get_vertices
{
    my $ritems = shift;

    my @array = @$ritems;
    my @new = ();

    shift(@array);
    shift(@array);
    shift(@array);
    my $nvert = shift(@array);

    foreach my $vert (@array) {
	push(@new,$vert);
	$nvert--;
	last unless $nvert;
    }

    die "*** Not enough vertices... @$ritems\n" if $nvert;

    return \@new;
}

###################################

sub delete_unused
{
    my ($self) = @_;

    $self->make_used();

    my $nodes = $self->{nodes};
    foreach my $node (values %$nodes) {
      unless( $node->{used} ) {
        $self->delete_node($node);
      }
    }
}

sub make_used
{
    my ($self) = @_;

    my $nodes = $self->{nodes};
    foreach my $node (values %$nodes) {
      $node->{used} = 0;
    }

    my $elems = $self->{elems};
    foreach my $item (values %$elems) {
      $self->inc_used($item->{vert});
    }

    my $lines = $self->{lines};
    foreach my $item (values %$lines) {
      $self->inc_used($item->{vert});
    }
}

sub inc_used
{
    my ($self,$verts) = @_;
    my $nodes = $self->{nodes};

    foreach my $number (@$verts) {
      my $node = $nodes->{$number};
      #print STDERR 
      $node->{used}++;
    }
}

###################################

sub set_verbose
{
    my ($self,$verbose) = @_;

    $self->{verbose} = $verbose;
}

###################################

sub test_grd
{
    my @files = @_;

    my $grid = new grd;

    print "reading...\n";
    $grid->readgrd($files[0]);
    print "writing...\n";
    $grid->writegrd;
}

###################################
#&test_grd(@ARGV);
###################################
1;
###################################

#!/usr/bin/perl -w
#
##############################################################
#
# version 1.3.1
#
# 26.08.2005            print_info, bug in make_unique
# 20.10.2005            version control introduced
# 24.10.2005            make_dxy
#
##############################################################
 
use strict;
 
package grdline;
 
##############################################################
 
sub new
{
    my $self;
 
    $self =     {
                         n           =>      0
                        ,x           =>      []
                        ,y           =>      []
                        ,closed      =>      -1
                        ,convex      =>      -1
                        ,xmin        =>      -1
                        ,xmax        =>      -1
                        ,ymin        =>      -1
                        ,ymax        =>      -1
			,dxy         =>      []
                };
 
    bless $self;
    return $self;
}
 
##############################################################

sub set_line {

  my ($self,$x,$y) = @_;

  my $nx = @$x;
  my $ny = @$y;

  die "Coordinates of different length: $nx $ny\n" if $nx != $ny;

  $self->{n} = $nx;

  my @xnew = @$x;
  my @ynew = @$y;

  $self->{x} = \@xnew;
  $self->{y} = \@ynew;

  $self->make_unique();
  $self->is_closed();
  $self->is_convex();
  $self->set_xy_min_max();
}

sub print_info {

  my ($self,$text,$verbose) = @_;

  print STDERR "$text\n";
  print STDERR "n: $self->{n}   closed: $self->{closed}" .
    "   convex: $self->{convex}\n";
  print STDERR "xy-min/max: $self->{xmin} $self->{ymin}" .
    " $self->{xmax} $self->{ymax}\n";
}

sub is_convex {

  my ($self) = @_;

  return $self->{convex} if $self->{convex} != -1;

  $self->{convex} = 0;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xl,$yl);
  my ($xm,$ym) = ($$x[$n-1],$$y[$n-1]);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xl,$yl) = ($xm,$ym);
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    return 0 unless( lefton($xl,$yl,$xm,$ym,$xn,$yn) );
  }

  $self->{convex} = 1;
  return 1;
}

sub in_convex {

  my ($self,$x0,$y0) = @_;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    #print STDERR "in_convex: $n $i $xm,$ym,$xn,$yn,$x0,$y0\n";
    return 0 unless( lefton($xm,$ym,$xn,$yn,$x0,$y0) );
  }

  return 1;
}

sub make_closed {
  my ($self) = @_;
  $self->{closed} = 1;
}

sub is_closed {

  my ($self) = @_;

  return $self->{closed} if $self->{closed} != -1;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  if( $$x[0] != $$x[$n] or $$y[0] != $$y[$n] ) {
    $self->{closed} = 0;
    return 0;
  }

  $self->{n} = $n;
  pop(@$x);
  pop(@$y);

  $self->{closed} = 1;
  return 1;
}

sub make_dxy {

  my ($self) = @_;

  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};
  my $dxy = $self->{dxy};

  for(my $i=0;$i<$n;$i++) {
    my $j = ($i+1) % $n;		#if not closed ignore last entry
    my $dx = $x->[$i] - $x->[$j];
    my $dy = $y->[$i] - $y->[$j];
    $dxy->[$i] = sqrt($dx*$dx+$dy+$dy);
  }
}

sub make_unique {

  my ($self) = @_;

  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};

  my @newx = ($x->[0]);
  my @newy = ($y->[0]);

  for( my $i=1; $i<$n; $i++) {
    if( $$x[$i] != $$x[$i-1] or $$y[$i] != $$y[$i-1] ) {
      push(@newx,$$x[$i]);
      push(@newy,$$y[$i]);
    } else {
      print STDERR "*** eliminating double point...\n";
    }
  }

  $self->{n} = @newx;
  $self->{x} = \@newx;
  $self->{y} = \@newy;
}

sub area {
 
  my $self = shift;

  my $n = $self->{n}-1;
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = $self->get_center_point();
  my $area = 0.;
 
  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    $area += areat($xm,$ym,$xn,$yn,$xc,$yc);
  }

  return $area;
}

sub get_center_point {
 
  my $self = shift;
  my $n = $self->{n};
  my $x = $self->{x};
  my $y = $self->{y};

  my ($xc,$yc) = (0,0);
 
  for( my $i=0;$i<$n;$i++ ) {
    $xc += $$x[$i];
    $yc += $$y[$i];
  }
 
  return ($xc/$n,$yc/$n);
}

#--------------------------------------------------------------------

sub smooth_line {

  my ($self,$sigma) = @_;

  $self->make_dxy() unless $self->{dxy};

  my $x = smooth_val($sigma,$self->{dxy},$self->{x},$self->{closed});
  my $y = smooth_val($sigma,$self->{dxy},$self->{y},$self->{closed});

  return ($x,$y);
}

sub smooth_val {

  my ($sigma,$dxy,$val,$is_closed) = @_;

  my ($imin,$imax,$j);
  my $pi = atan2(1,1);
  my $eps = 1.e-7;

  my $n = scalar @$dxy;
  my $a = 1/(sqrt(2*$pi)*$sigma);
  my $b = -1/(2*$sigma*$sigma);

  my @valnew = @$val;	#copy old values to new ones

  if( $is_closed ) {
    $imin = $n;
    $imax = 2*$n;
  } else {
    $imin = 1;
    $imax = $n-1;
  }

  for( my $i=$imin;$i<$imax;$i++) {
    my $tw = 0.;
    my $newval = 0.;
    my $jmin = $i - $n;
    my $jmax = $i + $n;
    unless( $is_closed ) {
      $jmin = 0; $jmax = $n-1;
    }

    my $kern_old = $a;
    my $tdist = 0.;
    $j = $i;
    while( $kern_old > $eps and $j++ < $jmax ) {
      my ($jt,$jo) = ($j%$n,($j-1)%$n);
      my $dist = $dxy->[$jo];
      $tdist += $dist;
      my $kern = $a * exp($b*$tdist*$tdist);
      my $w = 0.5 * $dist * ($kern_old+$kern);
      $tw += $w;
      $newval += 0.5 * $w * ($val->[$jt]+$val->[$jo]);
      $kern_old = $kern;
    }

    $kern_old = $a;
    $tdist = 0.;
    $j = $i;
    while( $kern_old > $eps and $j-- > $jmin ) {
      my ($jt,$jo) = ($j%$n,($j+1)%$n);
      my $dist = $dxy->[$jt];
      $tdist += $dist;
      my $kern = $a * exp($b*$tdist*$tdist);
      my $w = 0.5 * $dist * ($kern_old+$kern);
      $tw += $w;
      $newval += 0.5 * $w * ($val->[$jt]+$val->[$jo]);
      $kern_old = $kern;
    }

    $j = $i % $n;
    if( $tw > 0. ) {
      $valnew[$j] = $newval / $tw;
    }
  }
}

#--------------------------------------------------------------------

sub areat {

  my ($x1,$y1,$x2,$y2,$x3,$y3) = @_;

  return 0.5 * ( ($x2-$x1) * ($y3-$y1) - ($x3-$x1) * ($y2-$y1) );
}

sub lefton  { return ( areat(@_) >= 0 ); }
sub righton { return ( areat(@_) <= 0 ); }

sub angle {

  my ($x0,$y0,$xa,$ya,$xb,$yb) = @_;

  # angle between a-0-b (angle is on 0)

  my $dxa = $xa - $x0;
  my $dya = $ya - $y0;
  my $dxb = $xb - $x0;
  my $dyb = $yb - $y0;

  my $moda = sqrt( $dxa*$dxa + $dya*$dya );
  my $modb = sqrt( $dxb*$dxb + $dyb*$dyb );
  my $mod = $moda * $modb;

  if( $mod <= 0. ) {
    return 0.;
  }

  my $cos = (  $dxa * $dxb + $dya * $dyb ) / $mod;
  my $sin = ( -$dya * $dxb + $dxa * $dyb ) / $mod;

  if( $cos > 1. ) {
    $cos = 1.;
  }

  my $alpha = atan2( sqrt(1-$cos*$cos) , $cos );
  $alpha = -$alpha if( $sin < 0. );

  return $alpha;
}

sub get_min_max {
 
  my $a = shift;
 
  my ($min,$max);
 
  $min = $max = $$a[0];
  foreach my $c (@$a) {
    $max = $c if $c > $max;
    $min = $c if $c < $min;
  }
 
  return ($min,$max);
}

#--------------------------------------------------------------------

sub set_xy_min_max {
 
  my ($self) = @_;
 
  ($self->{xmin},$self->{xmax}) = get_min_max($self->{x});
  ($self->{ymin},$self->{ymax}) = get_min_max($self->{y});
}

sub in_min_max {

  my ($self,$x0,$y0) = @_;

  if( $self->{xmin} <= $x0 and $x0 <= $self->{xmax} ) {
    if( $self->{ymin} <= $y0 and $y0 <= $self->{ymax} ) {
	return 1;
    }
  }
  return 0;
}

sub in_line {

  my ($self,$x0,$y0) = @_;

  #print STDERR "in_line: $self->{n} $self->{convex} $self->{closed}\n";

  if( $self->{convex} ) {
    return $self->in_convex($x0,$y0);
  }

  return 0 unless $self->in_min_max($x0,$y0);

  return $self->in_any_line($x0,$y0);
}

sub in_any_line {

  my ($self,$x0,$y0) = @_;

  my $n = $self->{n} - 1;
  my $x = $self->{x};
  my $y = $self->{y};

  my $angle;
  my $pi = 4*atan2(1,1);

  my ($xm,$ym);
  my ($xn,$yn) = ($$x[$n],$$y[$n]);

  for (my $i=0;$i<=$n;$i++) {
    ($xm,$ym) = ($xn,$yn);
    ($xn,$yn) = ($$x[$i],$$y[$i]);

    $angle += angle($x0,$y0,$xm,$ym,$xn,$yn);
  }

  if( abs($angle) < $pi ) {
    return 0;
  } else {
    return 1;
  }
}

###################################
1;
###################################

