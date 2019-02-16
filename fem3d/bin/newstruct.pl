#!/usr/bin/perl -ws

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use lib ("$ENV{SHYFEMDIR}/fem3d/bin");

use strict;
use fortran;
 
my $fortran = new fortran;

my @ignore_routines = (
		         "subpar3.f"
			,"subnls.f"
			,"subdts.f"
		);

my @ignore_files = (
			"subomp_dummy.f"

			,"subgotm_mud.f"
			,"submud_dummy.f"

			,"subnetcdf_admin.f"
			,"subnetcdf_dummy.f"

			,"simsys_lp.f"
			,"simsys_pard.f"
			,"simsys_spk.f"

			,"newnohydro.f"
			,"nonhydro_lp.f"
			,"nonhydro_pard.f"
			,"nonhydro_spk.f"

			,"aquabc_fem_interface.f"
			,"ecological_dummy.f"
			,"subbfm.f"
			,"weutro.f"

			,"subnos.f"
			,"subous.f"
			,"subets.f"
			,"subsss.f"
		        ,"subpar3.f"
			,"subnls.f"
			,"subdts.f"
			,"subfil.f"
			,"nosutil.f"
			,"ousutil.f"
			,"genutil.f"

			,"submod.f"
			,"subcst.f"
		);

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

$fortran->check_uniqueness($::more) if $::unique;
$fortran->show_functions() if $::functions;
$fortran->show_routines() if $::routines;
$fortran->show_files() if $::files;

$fortran->show_implicit() if $::implicit;
$fortran->show_content() if $::content;

$fortran->calling_sequence($::struct) if $::struct;
$fortran->print_calling_times($::min_count);

$fortran->write_files($::quiet,$::all) if $::write;

#----------------------------------------------------------------

if( $::stats ) {

print "\n===================================================\n";
print "Statistics of routines found:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %{$fortran->{defined}}) ) {

  #print "........ $name\n";

  my $def = $fortran->{defined}->{$name};
  my $called = $fortran->{called}->{$name};
  my $calling = $fortran->{calling}->{$name};
  my $calledby = $fortran->{calledby}->{$name};
  my $item = $fortran->{routines}->{$name};
  my $type = $item->{type};
  my $file = $item->{file};

  $calling = defined($calling)  ? $calling  : "";
  $calling =~ s/^\s+//;
  my @f = split(/\s+/,$calling);
  my $nc = @f;

  $def = defined($def) ? $def : 0;
  $called = defined($called)  ? $called  : 0;
  print "$name  type: $type  defined: $def  calling: $nc  ";
  print "called: $called  file: $file\n";

  next unless $::more;

  my $list1 = $fortran->{calling}->{$name};
  $list1 = "" unless $list1;
  print "   is calling:   $list1\n";
  my $list2 = $fortran->{calledby}->{$name};
  $list2 = "" unless $list2;
  print "   is called by: $list2\n";
}

}

if( $::not_defined ) {

print "\n===================================================\n";
print "Programs called but not defined:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %{$fortran->{called}}) ) {
  my $def = exists($fortran->{defined}->{$name}) 
		? $fortran->{defined}->{$name} : 0;
  my $cal = exists($fortran->{called}->{$name})  
		? $fortran->{called}->{$name}  : 0;
  print "defined: $def    called: $cal    $name\n" if $def == 0;
}

}

if( $::not_called ) {

print "\n===================================================\n";
print "Programs defined but not called:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %{$fortran->{defined}}) ) {
  my $def = exists($fortran->{defined}->{$name}) 
		? $fortran->{defined}->{$name} : 0;
  my $cal = exists($fortran->{called}->{$name})  
		? $fortran->{called}->{$name}  : 0;
  print "defined: $def    called: $cal    $name\n" if $cal == 0;
}

}

if( $::compile ) {

print "\n===================================================\n";
print "Compiling programs:";
print "\n===================================================\n\n";

make_dummy();

my $file_list = {};

my @flist = $fortran->get_file_list();
foreach my $file (@flist) {
  my $fitem = $fortran->get_file_item($file);
  my $nlines = $fitem->{nlines};
  my $has_program = $fitem->{has_program};
  #print STDERR "$file  $nlines  $has_program\n";
  next if $has_program;
  make_compile($fortran,$file,$file_list);
}

print "\n===================================================\n";
print "Compiling dependencies:";
print "\n===================================================\n\n";

foreach my $file (sort keys %$file_list) {
  my $flist = $file_list->{$file};
  print "$file: $flist\n";
}

}

#--------------------------------------------------------------

sub make_dummy {

  open(FILE,">main_dummy.f");
  print FILE "\tend\n";
  close(FILE);
}

sub make_compile {

  my ($fortran,$file,$file_list) = @_;

  my @output = `gfortran main_dummy.f $file 2>&1`;
  #print "=======================================\n";
  #print "compilation output for file $file\n";
  #print "=======================================\n";
  #print @output;

  my ($r,$f) = parse_compile_output($fortran,$file,\@output);

  print "compiling $file:\n";
  my $rlist = join(" ",@$r);
  print "  routines: $rlist\n";
  my $flist = join(" ",@$f);
  print "  files: $flist\n";

  $file_list->{$file} = $flist;
}

sub parse_compile_output {

  my ($fortran,$file,$output) = @_;

  my %routines = ();
  my %files = ();

  foreach my $line (@$output) {
    if( $line =~ /undefined reference to \`(\w+)_\'/ ) {
      my $routine = $1;
      $routines{$routine}++;
      my $ritem = $fortran->get_routine_item($routine);
      my $file = $ritem->{file};
      $file = "no_info" unless $file;
      #print "^^^^^^^^^^^^ $routine: $file  $type  $name\n";
      $files{$file}++;
    }
  }

  my @r = sort keys %routines;
  my @f = sort keys %files;

  return (\@r,\@f);
}

#--------------------------------------------------------------

