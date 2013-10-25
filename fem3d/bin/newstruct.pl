#!/usr/bin/perl -w

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use lib ("$ENV{SHYFEMDIR}/fem3d/bin");

use strict;
use fortran;
 
my $fortran = new fortran;

my @ignore_files = (
		         "subpar.f"
			,"subfnm.f"

			,"subsse.f"

			,"subomp_dummy.f"

			,"subgotm_mud.f"
			,"submud_dummy.f"

			,"subnetcdf_admin.f"
			,"subnetcdf_dummy.f"

			,"simsys_lp.f"
			,"simsys_pard.f"
			,"simsys_spk.f"

			,"nonhydro_lp.f"
			,"nonhydro_pard.f"
			,"nonhydro_spk.f"

			,"aquabc_fem_interface.f"
			,"ecological_dummy.f"
			,"subbfm.f"
			,"weutro.f"
		);

my $show_stats = 0;
my $not_defined = 0;
my $not_called = 0;
my $quiet = 0;

$fortran->{no_sub_in_program} = 1;	# no insert of sub in files with main

$fortran->set_ignore_files(\@ignore_files);
$fortran->parse_files($quiet);
$fortran->check_uniquness();
$fortran->parse_routines();

#$fortran->show_functions();
$fortran->show_routines();

if( $show_stats ) {

print "\n===================================================\n";
print "Statistics of routines found:";
print "\n===================================================\n\n";

foreach my $name ( sort (keys %{$fortran->{defined}}) ) {

  print "........ $name\n";

  my $def = $fortran->{defined}->{$name};
  my $called = $fortran->{called}->{$name};
  my $calling = $fortran->{calling}->{$name};
  my $calledby = $fortran->{calledby}->{$name};
  my $item = $fortran->{routines}->{$name};
  my $type = $item->{type};

  $def = defined($def) ? $def : 0;
  $called = defined($called)  ? $called  : 0;
  print "$name   defined: $def   called: $called  type: $type\n";
  my $list1 = $fortran->{calling}->{$name};
  $list1 = "" unless $list1;
  print "   is calling:   $list1\n";
  my $list2 = $fortran->{calledby}->{$name};
  $list2 = "" unless $list2;
  print "   is called by: $list2\n";
}

}

if( $not_defined ) {

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

if( $not_called ) {

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

