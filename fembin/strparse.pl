#!/usr/bin/perl -s -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# parses STR file and writes info to stdout or grd file
#
# possible command line options: see subroutine FullUsage
#
# version       1.7     02.11.2023     -restart implemented
# version       1.6     11.06.2023     -essential implemented
# version       1.5     09.06.2023     handle flux/extra section with descrip
# version       1.4     07.06.2023     use ScriptPath() to determine path
# version       1.3     09.05.2023     handle also boxes and only bas files
# version       1.2     30.03.2023     restructured and -collect, -reduce
# version       1.1     24.03.2023     handle extra and flux sections
# version       1.0     ?              old version, partially functional
#
#--------------------------------------------------------

use FindBin qw($Bin); use lib "$Bin/../lib/perl"; 

use str;
use grd;
use strict;

#-------------------------------------------------------------
# command line options
#-------------------------------------------------------------
$::h = 0 unless $::h;
$::help = 0 unless $::help;
$::quiet = 0 unless $::quiet;
$::bnd = 0 unless $::bnd;
$::files = 0 unless $::files;
$::extra = 0 unless $::extra;
$::flux = 0 unless $::flux;
$::zip = 0 unless $::zip;
$::rewrite = 0 unless $::rewrite;
$::restart = 0 unless $::restart;
$::essential = 0 unless $::essential;
$::sect = "" unless $::sect;
$::txt = 0 unless $::txt;
$::debug = 0 unless $::debug;
$::value = "" unless $::value;
$::replace = "" unless $::replace;
$::simtime = 0 unless $::simtime;
$::collect = "" unless $::collect;
$::reduce = "" unless $::reduce;
#-------------------------------------------------------------

#-------------------------------------------------------------
# info on names in sections
#-------------------------------------------------------------
@::bound_names = qw/ boundn conzn tempn saltn vel3dn
		bio2dn sed2dn tox3dn
		bfm1bc bfm2bc bfm3bc /;
@::name_names = qw/ bound wind rain qflux ice restrt gotmpa
		bio bios toxi
		conzin saltin tempin zinit /;
@::aquabc_names = qw/ biocon bioscon biolight bioaow 
		bioaos bioph biotemp bioload /;
@::lagrg_names = qw/ lagra lagrf lagrt /;
@::sedtr_names = qw/ sedp sedt sedcon /;
#-------------------------------------------------------------

my $file = $ARGV[0];

$::scriptpath = ScriptPath();

my $str = new str;
$str->{quiet} = $::quiet;
$str->{nocomment} = $::essential;
$str->read_str($file) if $file;

if( $::h or $::help ) {
  FullUsage();
} elsif( not $file ) {
  Usage();
} elsif( $::bnd ) {
  show_bnd_nodes($str);
} elsif( $::extra ) {
  show_extra_nodes($str);
} elsif( $::flux ) {
  show_flux_nodes($str);
} elsif( $::files ) {
  show_files($str);
} elsif( $::zip ) {
  my $files = show_files($str);
  push(@$files,$file);		#add str-file to archive
  zip_files($files);
} elsif( $::rewrite ) {
  #$str->print_sections();;
  $str->write_str("new.str");;
  print STDERR "str-file written to new.str\n";
} elsif( $::restart ) {
  prepare_restart($str);
  $str->write_str("restart.str");;
  print STDERR "str-file written to restart.str\n";
} elsif( $::sect ) {
  show_sect($str,$::sect);
} elsif( $::simtime ) {
  my $val = "";
  $val = show_value($str,"itanf");
  print "$val\n";
  $val = show_value($str,"itend");
  print "$val\n";
} elsif( $::value ) {
  if( $::replace ) {
    replace_value($str,$::value,$::replace);
    $str->write_str("replace.str");;
    print STDERR "str-file written to replace.str\n";
  } else {
    my $val = show_value($str,$::value);
    $val = "(unknown)" unless $val;
    #print "$::value = $val\n";
    print "$val\n";
  }
} elsif( $::collect ) {
  if( $::collect eq "1" ) {
    print STDERR "please specify directory for collecting: -collect=dir\n";
    exit 1;
  }
  print STDERR "copy simulation of str file into directory $::collect\n";
  collect_str($str,$::collect);
} else {
  Usage();
}

#------------------------------------------------------------

sub Usage {

  print STDERR "Usage: strparse.pl [-h|-help] [-options] str-file\n";
  exit 0;
}

sub FullUsage {

  print STDERR "Usage: strparse.pl [-h|-help] [-options] str-file\n";
  print STDERR "  options:\n";
  print STDERR "    -h!-help      this help screen\n";
  print STDERR "    -quiet        be as quiet as possible\n";
  print STDERR "    -bnd          extract boundary nodes\n";
  print STDERR "    -files        extract names of forcing files\n";
  print STDERR "    -extra        extract extra nodes\n";
  print STDERR "    -flux         extract flux nodes\n";
  print STDERR "    -zip          zips forcing files, grid, str in one file\n";
  print STDERR "    -rewrite      rewrite the str file\n";
  print STDERR "    -restart      create a restart from the str file\n";
  print STDERR "    -essential    only write active sections, no comments\n";
  print STDERR "    -value=var    show value of var ([sect:]var)\n";
  print STDERR "    -replace=val  replace value of var with val and rewrite\n";
  print STDERR "    -simtime      prints start and end time of simulation\n";
  print STDERR "    -collect=dir  copies str with data into dir\n";
  print STDERR "    -reduce       reduce data files for collect, else copy\n";
  #print STDERR "    -sect=sect    writes contents of section\n";
  print STDERR "    -txt          write nodes as text and not in grd format\n";
  print STDERR "  if -replace is given also -value must be specified\n";
  print STDERR "  using -collect=dir -reduce choses minimal data for sim\n";
  exit 0;
}

#------------------------------------------------------------

sub ScriptPath {

  my $script = `realpath $0`;
  chomp $script;
  my $scriptpath = `dirname $script`;
  chomp $scriptpath;

  #print STDERR "$script\n";
  #print STDERR "$scriptpath\n";

  return $scriptpath;
}

sub show_flux_nodes {

  my $str = shift;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  my $basin = $str->get_basin();
  my $grid = new grd;
  $grid->readgrd("$basin.grd");

  unless( $str->has_section("flux") ) {
    die "*** no flux section found...\n";
  }

  open_nodes_file("flux_str");

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "flux" ) {
      my ($rlist,$dlist) = get_node_list($str,$sect);
      my $slist = make_sections($rlist);
      my $n = 0;
      foreach my $ra (@$slist) {
        $n++;
        write_line($n,-1,$grid,$ra);
	my $descr = shift(@$dlist);
        print_section($n,$ra,$descr);
      }
    }
  }

  close_nodes_file();
}

sub print_section {

  my ($n,$ra,$descr) = @_;

  $descr = "" unless $descr;
  my $na = @$ra;

  print STDERR "section $n ($na)  $descr\n";
  foreach my $n (@$ra) {
    print STDERR "  $n";
  }
  print STDERR "\n";
}

sub make_sections {

  my $rlist = shift;

  my @slist = ();
  my @list = ();

  foreach my $n (@$rlist) {
    if( $n == 0 ) {
      my $nl = @list;
      next if $nl == 0;
      my @new = @list;
      push(@slist,\@new);
      @list = ();
    } else {
      push(@list,$n);
    }
  }

  return \@slist;
}

sub show_extra_nodes {

  my $str = shift;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  my $basin = $str->get_basin();
  my $grid = new grd;
  $grid->readgrd("$basin.grd");

  unless( $str->has_section("extra") ) {
    die "*** no extra section found...\n";
  }

  open_nodes_file("extra_str");

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "extra" ) {
      my $ra = $sect->{array};
      my ($rlist,$dlist) = get_node_list($str,$sect);
      write_nodes(-1,$grid,$rlist);
      print_nodes($rlist,$dlist);
    }
  }

  close_nodes_file();
}

sub get_node_list {

  my ($str,$sect) = @_;

  my $sect_name = $sect->{name}; 

  return($sect->{array},$sect->{description});
}

sub print_nodes {

  my ($rlist,$dlist) = @_;

  my $na = @$rlist;
  my $nd = @$dlist;

  print STDERR "node list: $na $nd\n";

  foreach my $n (@$rlist) {
    print STDERR "$n ";
  }
  print STDERR "\n";

  return unless $nd;

  foreach my $n (@$rlist) {
    my $descr = shift(@$dlist);
    print STDERR "  $n   $descr\n";
  }
}

sub show_bnd_nodes {

  my $str = shift;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  my $basin = $str->get_basin();
  my $grid = new grd;
  $grid->readgrd("$basin.grd");

  open_nodes_file("bnd_str");

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    if( $sect->{name} eq "bound" ) {
      my ($itype,$rlist) = parse_bnd_nodes($str,$sect);
      write_nodes($itype,$grid,$rlist);
    }
  }

  close_nodes_file();
}

sub parse_bnd_nodes {

  my ($str,$sect) = @_;

  my $sect_name = $sect->{name}; 
  my $sect_number = $sect->{number}; 
  my $sect_id = $sect->{id}; 

  my @list = ();

  my $ibtyp = $str->get_value('ibtyp',$sect_name,$sect_number);
  $ibtyp = 1 if not defined $ibtyp;
  my $value = $str->get_value('kbound',$sect_name,$sect_number);
  if( defined $value ) {
    if( ref($value) eq "ARRAY" ) {
      write_array_new($value,"$sect_id :  kbound = \n");
      @list = @$value;
    } else {
      print "$sect_id :    kbound = $value\n";
      push(@list,$value);
    }
  }

  return($sect_number,\@list);
}

#------------------------------------------------------------

sub open_nodes_file {

  my $name = shift;

  if( $::txt ) {
    $::outfile = "$name.txt";
  } else {
    $::outfile = "$name.grd";
  }
  open(OUT,">$::outfile");
}

sub close_nodes_file {

  close(OUT);
  print STDERR "nodes written to file $::outfile\n";
}

sub write_line {

  my ($nline,$itype,$grid,$rlist) = @_;

  if( $::txt ) {
    write_line_to_txt($itype,$rlist);
  } else {
    write_line_to_grd($nline,$itype,$grid,$rlist);
  }
}

sub write_line_to_txt {

  my ($id,$list) = @_;

  foreach my $val (@$list) {
    print OUT "$val ";
  }
  print OUT "\n";
}

sub write_line_to_grd {

  my ($nline,$ibtyp,$grid,$list) = @_;

  my $nval = 5;                 # how many values on one line

  my $n = 0;
  #return if( $ibtyp == 0 );
  my $type = $ibtyp;

  foreach my $node (@$list) {
    my $item = $grid->get_node($node);
    my $x = $item->{x};
    my $y = $item->{y};
    $n++;
    $type = $n if $ibtyp < 0;
    print OUT "1 $node 0 $x $y\n";
  }
  my $j = 0;
  $type = $n if $ibtyp < 0;
  print OUT "3 $nline $type $n\n";
  foreach my $nn (@$list) {
    $j++;
    print OUT "  $nn";
    print OUT "\n" if $j%$nval == 0;
  }
  print OUT "\n" if not $j%$nval == 0;
}

sub write_nodes {

  my ($itype,$grid,$rlist) = @_;

  if( $::txt ) {
    write_nodes_to_txt($itype,$rlist);
  } else {
    write_nodes_to_grd($itype,$grid,$rlist);
  }
}

sub write_nodes_to_txt {

  my ($id,$list) = @_;

  foreach my $val (@$list) {
    print OUT "$val\n";
  }
}

sub write_nodes_to_grd {

  my ($ibtyp,$grid,$list) = @_;

  my $n = 0;
  #return if( $ibtyp == 0 );
  my $type = $ibtyp;

  foreach my $node (@$list) {
    my $item = $grid->get_node($node);
    my $x = $item->{x};
    my $y = $item->{y};
    $n++;
    $type = $n if $ibtyp < 0;
    print OUT "1 $node $type $x $y\n";
  }
}

#------------------------------------------------------------

sub prepare_restart {

  my $str = shift;

  my $simul = $str->get_simul();
  replace_value($str,'itrst',-1,'para');
  replace_value($str,'restrt',"'${simul}.rst'",'name');
  $simul .= "_rst";
  $str->set_simul($simul);
}

sub show_files {

  my $str = shift;

  my @files = ();
  my @items = ();
  my $newfile = "";
  my $newitem = "";

  my $basin = $str->get_basin();
  $basin =~ s/^\s+//;
  if( file_exists("$basin.grd") ) {
    $basin .= ".grd";
  } elsif( file_exists("$basin.bas") ) {
    $basin .= ".bas";
  } else {
    die "*** cannot find basin $basin: no .grd or .bas file found";
  }
  print "  basin :    grid = '$basin'\n";
  push(@files,$basin);

  my %item = ();
  $item{"section"} = "title";
  $item{"varname"} = "basin";
  $item{"value"} = $basin;
  push(@items,\%item);

  if( file_exists("boxes.txt") ) {
    my %bitem = ();
    print "  boxes :    boxes = 'boxes.txt'\n";
    push(@files,"boxes.txt");
    $bitem{"section"} = "none";
    $bitem{"varname"} = "boxes";
    $bitem{"value"} = "boxes.txt";
    push(@items,\%bitem);
  }

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};
    my $section_id = $sect->{id};

    $newfile = "";
    $newitem = "";
    if( $sect->{name} eq "bound" ) {
      ($newfile,$newitem) = show_name($str,$sect,\@::bound_names);
    } elsif( $sect->{name} eq "name" ) {
      ($newfile,$newitem) = show_name($str,$sect,\@::name_names);
    } elsif( $sect->{name} eq "aquabc" ) {
      ($newfile,$newitem) = show_name($str,$sect,\@::aquabc_names);
    } elsif( $sect->{name} eq "lagrg" ) {
      ($newfile,$newitem) = show_name($str,$sect,\@::lagrg_names);
    } elsif( $sect->{name} eq "sedtr" ) {
      ($newfile,$newitem) = show_name($str,$sect,\@::sedtr_names);
    }

    push(@files,@$newfile) if $newfile;
    push(@items,@$newitem) if $newitem;
  }

  return (\@files,\@items);
}

sub replace_value {

  my ($str,$name,$replace,$sect_name,$sect_number) = @_;

  #$sect_name = "para" unless $sect_name;
  #$sect_number = "" unless $sect_number;

  $str->set_value($name,$replace,$sect_name,$sect_number);
}

sub show_value {

  my ($str,$name,$sect_name,$sect_number) = @_;

  #$sect_name = "para" unless $sect_name;
  #$sect_number = "" unless $sect_number;

  return $str->get_value($name,$sect_name,$sect_number);
}

sub show_name {

  my ($str,$sect,$var_names) = @_;

  my @files = ();
  my @items = ();

  my $sect_name = $sect->{name}; 
  my $sect_number = $sect->{number}; 
  my $sect_id = $sect->{id}; 

  if( $::debug ) {
    print STDERR "searching files: $sect_name ($sect_number) $sect_id\n";
  }

  foreach my $name (@$var_names) {
    my $value = $str->get_value($name,$sect_name,$sect_number);
    print STDERR "  $name  (no value)\n" if not defined $value and $::debug;
    if( defined $value ) {
      print STDERR "  $name  $value\n" if $::debug;
      print "  $sect_id :    $name = $value\n";
      $value =~ s/\'//g;
      push(@files,$value);
      my %item = ();
      $item{"section"} = $sect_id;
      $item{"varname"} = $name;
      $item{"value"} = $value;
      push(@items,\%item);
    }
  }

  return (\@files,\@items);
}

sub write_array_new {

  my ($array,$extra) = @_;

  my $nval = 10;		# how many values on one line

  print "$extra" if $extra;

  my $i = 0;
  foreach my $item (@$array) {
    $i++;
    print " $item";
    print "\n" if $i%$nval == 0;
  }
  print "\n" unless $i%$nval == 0;
}

sub write_txt {

  my ($id,$value) = @_;

  if( ref($value) ne "ARRAY" ) {
    my @value = ();
    $value[0] = $value;
    $value = \@value;
  }

  my $n = @$value;
  print "$id  $n\n";
  while( $n-- ) {
    my $val = shift(@$value);
    print "$val\n";
  }
}

sub show_sect {

  my ($str,$sname) = @_;

  my $sections = $str->{sections};
  my $sequence = $str->{sequence};

  foreach my $section (@$sequence) {
    my $sect = $sections->{$section};

    my $nnn = $sect->{name}; print "section: $nnn\n";

    if( $sect->{name} eq "$sname" ) {
        #show_name($str,$sect,\@::bound_names);
	my $array = $sect->{array};
	foreach my $numb (@$array) {
	  print "$numb\n";
	}
    }
  }
}

sub zip_files {

  my ($files) = @_;

  my $zipfile = "new_str_zip.zip";

  unlink "$zipfile" if ( -f $zipfile );

  foreach my $file (@$files) {
    print STDERR "zip file: $file\n";
    system("zip $zipfile $file");
  }

  print STDERR "files have been zipped into file $zipfile\n";
}

sub collect_str {

  my ($str,$dir) = @_;

  my $basin;

  my $itanf = show_value($str,"itanf");
  my $itend = show_value($str,"itend");
  my $date = show_value($str,"date");
  print "  time interval: $itanf - $itend ($date)\n";

  print STDERR "handling the following files:\n";

  my ($files,$items) = show_files($str);

  my $newinput = "input";

  my $input = "$dir/$newinput";
  system("mkdir -p $input");

  print STDERR "copying/condensing the following files:\n";

  foreach my $item (@$items) {
    my $file = $item->{"value"};
    my $filename = $file;
    $filename =~ s/.*\///;
    if( $file =~ /\.grd$/ ) {
      print STDERR "  copying grd file $file to $dir\n";
      system("cp $file $dir");
      $item->{"newdir"} = "";
    } elsif( $file =~ /\.bas$/ ) {
      print STDERR "  copying bas file $file to $dir\n";
      system("cp $file $dir");
      $item->{"newdir"} = "";
    } elsif( $file =~ /\.nml$/ ) {
      print STDERR "  copying nml file $file to $dir\n";
      system("cp $file $dir");
      $item->{"newdir"} = "";
    } elsif( $file =~ /boxes.txt$/ ) {
      print STDERR "  copying boxes file $file to $dir\n";
      system("cp $file $dir");
      $item->{"newdir"} = "";
    } elsif( not $::reduce ) {
      print STDERR "  copying file $file to $input\n";
      system("cp $file $input");
      $item->{"newdir"} = "$newinput/";
    } else {
      my $type = `$::scriptpath/shyfile $file`;
      chomp($type);
      if( $type eq "FEM" ) {
        print STDERR "  reducing ($type) $file to $input as $filename\n";
        my $command = "femelab -tmin $itanf -tmax $itend";
        my $options = "-out -inclusive";
        system("[ -f out.fem ] && rm -f out.fem");
        system("$command $options $file > /dev/null");
        system("mv out.fem $input/$filename");
      } elsif( $type eq "TS" ) {
        print STDERR "  reducing ($type) $file to $input as $filename\n";
        my $command = "tselab -tmin $itanf -tmax $itend";
        my $options = "-out -inclusive";
        system("[ -f out.txt ] && rm -f out.fem");
        system("$command $options $file > /dev/null");
        system("mv out.txt $input/$filename");
      } else {
        print STDERR "  cannot reduce type $type $file ...copying\n";
        system("cp $file $input");
      }
      $item->{"newdir"} = "$newinput/";
    }
    $item->{"filename"} = $filename;
  }

  print STDERR "changing the following names in the STR file\n";

  foreach my $item (@$items) {
    my $file = $item->{"value"};
    my $filename = $item->{"filename"};
    my $section = $item->{"section"};
    my $varname = $item->{"varname"};
    my $newdir = $item->{"newdir"};

    print "  $section $varname $filename\n";

    if( $varname eq "basin" ) {
      $filename =~ s/\.grd//;
      $basin = $filename;
      $str->set_basin($filename);
    } elsif( $section eq "none" ) {
      #nothing
    } else {
      replace_value($str,$varname,"'$newdir$filename'",$section);
    }
  }

  print STDERR "creating Makefile\n";

  create_makefile($dir,$basin);
  
  $str->write_str("$dir/$dir.str");;

  print STDERR "collected files are in directory $dir with $dir.str\n";
  print STDERR "you can pack all files by running one of the following:\n";
  print STDERR "  zip -r $dir.zip $dir/*\n";
  print STDERR "  tar cvzf $dir.tar.gz $dir\n";
  print STDERR "you can run the simulation by doing the following:\n";
  print STDERR "  cd $dir; make basin; make run\n";
}

sub file_exists {

  my $file = shift;

  if( -e $file ) {
    return 1;
  } else {
    return 0;
  }
}

sub create_makefile {

  my ($dir,$basin) = @_;

  open(MAKE,">$dir/Makefile");

  print MAKE "\n";
  print MAKE "BASIN = $basin.grd\n";
  print MAKE "STR = $dir.str\n";

my $var = <<'EOF';

default:

basin:
	shypre $(BASIN)

run:
	shyfem $(STR)

clean:
	-rm -f errout.dat

cleanall: clean
	-rm -f *.bas
	-rm -f *.shy *.rst *.flx *.ext *.log *.inf
	-rm -f boxes_*.txt

EOF

  print MAKE "$var\n";
}

