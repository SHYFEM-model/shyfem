#!/usr/bin/perl -ws

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use lib ("$ENV{SHYFEMDIR}/fem3d/bin");

use strict;
use fortran;
 
my $fortran = new fortran;

my @ignore_routines = ();
my @ignore_files = ();

#----------------------------------------------------------------

$::quiet = 0 unless $::quiet;
$::unique = 0 unless $::unique;
$::functions = 0 unless $::functions;
$::routines = 0 unless $::routines;
$::files = 0 unless $::files;
$::write = 0 unless $::write;
$::stats = 0 unless $::stats;
$::more = 0 unless $::more;
$::all = 0 unless $::all;
$::not_defined = 0 unless $::not_defined;
$::not_called = 0 unless $::not_called;
$::struct = "" unless $::struct;
$::no_main = 0 unless $::no_main;
$::min_count = 0 unless $::min_count;
$::debug = 0 unless $::debug;
$::parse_later = 0 unless $::parse_later;
$::subst = 0 unless $::subst;
$::clean = 0 unless $::clean;

$::implicit = 0 unless $::implicit;
$::content = 0 unless $::content;
$::compile = 0 unless $::compile;

#----------------------------------------------------------------

$fortran->{no_main} = $::no_main;	# do not insert subs in files with main
$fortran->{implicit} = $::implicit;
$fortran->{debug} = $::debug;
$fortran->{parse_later} = $::parse_later;

$fortran->set_ignore_files(\@ignore_files);
$fortran->set_ignore_routines(\@ignore_routines);

$fortran->parse_files($::quiet);
$fortran->parse_routines() if $::parse_later;

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

if( $::subst ) {
  subst_konst($fortran);
} elsif( $::clean ) {
  clean_files($fortran);
}

sub subst_konst {
  my $fortran = shift;
  subst_common($fortran,"mkonst","mkonst.h");
  subst_common($fortran,"pkonst","pkonst.h");
  subst_common($fortran,"nkonst","nbasin.h");
}
sub subst_femtim {
  my $fortran = shift;
  subst_common($fortran,"femtim","femtime.h");
  subst_common($fortran,"femtimu","femtime.h");
}

#subst_common($fortran,"saltv","baroc.h");
#subst_common($fortran,"cnv","conz.h");

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

$fortran->write_files($::quiet,$::all) if $::write;

#----------------------------------------------------------------

sub subst_common {

  my ($fortran,$common,$include) = @_;

  print STDERR "  substituting common /$common/ with include $include\n";

  my @list = sort keys %{$fortran->{all_routines}};
  foreach my $rname (@list) {
    my $ritem = $fortran->{all_routines}->{$rname};
    my $name = $ritem->{name};
    my $file = $ritem->{file};

    if( $ritem->{common}->{$common} ) {
      print STDERR "    routine $name contains common... substituting\n";
      substitute_common($fortran,$ritem,$common,$include);
      $fortran->set_changed($file);
    }
  }
}

sub substitute_common {

  my ($fortran,$ritem,$common,$include) = @_;

  my ($name,$list);
  my $code = $ritem->{code};
  my @new = ();
  my $in_specification = 1;
  my $ncode = @$code;
  my $nline = 0;

  my $alist = $ritem->{common}->{$common};
  my $hlist = {};
  foreach my $var (@$alist) {
    $var =~ s/\(.+$//;
    $hlist->{$var}++;
  }
  my $clist = join(",",@$alist);
  print STDERR "    list to be treated: $clist\n";

  foreach my $l (@$code) {		# common must not be continuation line
    $nline++;
    if( $in_specification and $nline > 1 and not is_conti($l) ) {
      my $line = $fortran->clean_line($l);
      if( $line ) {
        $line =~ s/\s+//g;
        #print "++++++++++++++ $line\n";
        if( my $what = $fortran->is_specification($line,$ritem) ) {
          #print "============== $what\n";
	  if( $what eq "common" ) {
  	    my ($found,$new) = treat_common($fortran,$common,$line);
	    if( $found ) {
    	      push(@new,$new) if $new;
	      unless( $ritem->{include}->{$include} ) {
    	        push(@new,"\tinclude '$include' !COMMON_GGU_SUBST");
	      }
	      $l = "COMMON_GGU_DELETED$l";
	    }
	  } elsif( $what eq "save" ) {
  	    my ($found,$new) = treat_save($fortran,$common,$line);
	    if( $found ) {
    	      push(@new,$new) if $new;
	      $l = "COMMON_GGU_DELETED$l";
	    }
	  } elsif( $what =~ /^declaration-(\w+)/i ) {
	    my $type = $1;
  	    my ($found,$new) = treat_declaration($fortran,$common
				,$line,$type,$hlist);
	    if( $found ) {
    	      push(@new,$new) if $new;
	      $l = "COMMON_GGU_DELETED$l";
	    }
	  }
        } else {
          #print "is no spec_ $l\n";
          $in_specification = 0;
        }
      }
    }
    push(@new,$l);
  }

  $ritem->{code} = \@new;
}

sub is_conti {

  my $line = shift;

  if( $line =~ /^     \S/ ) {
    return 1;
  } else {
    return 0;
  }
}

sub treat_declaration {

  my ($fortran,$common,$line,$type,$hlist) = @_;

  $line =~ s/^$type//;
  my $new = "";
  my $found = 0;

  my $vars = $fortran->split_var($line);
  #print "      declaration: $type - $line\n";
  foreach my $var (@$vars) {
    my $orig = $var;
    $var =~ s/\(.+$//;
    if( $hlist->{$var} ) {
      $found = 1;
    } else {
      $new .= "$orig, "
    }
  }

  if( $found ) {
    print STDERR "      declaration found: $type $line\n";
  }
  if( $new ) {
    $new =~ s/,\s*$//;
    $new = "\t$type " . $new . " !COMMON_GGU_SUBST";
  }

  return ($found,$new);
}

sub treat_save {

  my ($fortran,$common,$line) = @_;

  $line =~ s/^save//;
  my $new = "";
  my $found = 0;

  my $vars = $fortran->split_var($line);
  foreach my $var (@$vars) {
    if( $var eq "/$common/" ) {
      print STDERR "      save found: $var\n";
      $found = 1;
    } else {
      $new .= "$var, "
    }
  }
  if( $new ) {
    $new =~ s/,\s*$//;
    $new = "\tsave " . $new . " !COMMON_GGU_SUBST";
  }

  return ($found,$new);
}

sub treat_common {

  my ($fortran,$common,$line) = @_;

  my ($name,$list);
  $line =~ s/^common//;
  my $new = "";
  my $found = 0;

  while(1) { 
    ($name,$list,$line) = $fortran->next_common($line);
    last unless $name;
    if( $name eq $common ) {
      print STDERR "      common found: /$name/ $list\n";
      $found = 1;
    } else {
      $new .= "/$name/$list, "
    }
  }
  if( $new ) {
    $new =~ s/,\s*$//;
    $new = "\tcommon " . $new . " !COMMON_GGU_SUBST";
  }

  return ($found,$new);
}

#--------------------------------------------------------------

sub clean_files {

  my ($fortran) = @_;

  foreach my $filename (keys %{$fortran->{files}}) {
    clean_file($fortran,$filename)
  }
}

sub clean_file {

  my ($fortran,$filename) = @_;

  my $fitem = $fortran->{files}->{$filename};
  my $sequence = $fitem->{sequence};

  my $changed = 0;
  my @new = ();
  foreach my $ritem (@$sequence) {
    my $code = $ritem->{code};
    foreach my $line (@$code) {
      if( $line =~ /^COMMON_GGU_DELETED/ ) { $changed++; next; }
      if( $line =~ s/\s*!COMMON_GGU_SUBST\s*$// ) { $changed++; }
      push(@new,$line);
    }
  }

  return unless $changed;

  my $newfile = $filename . ".clean";
  print STDERR "cleaning file $filename ($changed) to $newfile\n";
  open(NEW,">$newfile");
  foreach my $line (@new) {
    print NEW "$line\n";
  }
  close(NEW);
}

