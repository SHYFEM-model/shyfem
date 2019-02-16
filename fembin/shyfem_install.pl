#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# inserts new shyfem section into dot file
#
#-----------------------------------------------------

$reset="";
$femdir = shift;
if( $femdir eq "-reset" ) {
  $reset = "-reset";
  $femdir = shift;
}

$in_section = 0;

unless( $femdir ) {
  print STDERR "shyfem_install.pl should not be run directly... aborting\n";
  exit 1;
}

while(<>) {

  chomp;

  if( /^\#--- SHYFEM --- set up start\s*$/ ) {
    $in_section = 1;
  }

  if( $in_section != 1 ) {
    print "$_\n";
  }

  if( /^\#--- SHYFEM --- set up end\s*$/ ) {
    write_section();
    $in_section = 2;
  }
}

if( $in_section == 0 ) {	# no section found
  write_section();
} elsif( $in_section == 2 ) {	# section found and written
  ;
} else {
  print STDERR "error in changing dot file... aborting";
  exit 1;
}

exit 0;

#---------------------------------------------------------------

sub write_section {

  return if $reset;		# reset - do not write anything

  print "\#--- SHYFEM --- set up start\n";
  print "export SHYFEM_INSTALL=$femdir\n";
  print "export SHYFEMDIR=\$SHYFEM_INSTALL\n";
  print "if [ -f \$SHYFEMDIR/fembin/shyfem_profile.sh ]; then\n";
  print	"  . \$SHYFEMDIR/fembin/shyfem_profile.sh\n";
  print "fi\n";
  print "\#--- SHYFEM --- set up end\n";
}
  
