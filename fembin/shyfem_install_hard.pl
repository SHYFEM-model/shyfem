#!/usr/bin/perl -s
#
# installs shyfem through hard references
#
# called by shell script
#
#---------------------------------------------

$::verbose = 0;

$reset = 0 unless $reset;
$install = 0 unless $install;

$dir = shift;

if( $reset ) {
  ;
} elsif( $install ) {
  ;
} else {
  die "one of -reset or -install must be set... aborting\n";
}

$change = 0;

if( $reset ) {
  reset_fem();
} else {
  install_fem();
}

exit $change;

#-------------------------------------------------------------

sub install_fem {

 while(<>) {

  if( /SHYFEMDIR/ ) {
    if( /^FEMDIR=/ ) {
      $change = 1;
      print "SHYFEMDIR=$dir        # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (shell): $_" if $::verbose;
    } elsif( /^use lib/ and /femlib\/perl/ ) {
      $change = 1;
      print "#ggu_hard_install $_";		# write old first
      print "use lib (\"$dir/femlib/perl\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (femlib/perl): $_" if $::verbose;
      next;
    } elsif( /^use lib/ and /ENV/ and /fembin/ ) {
      $change = 1;
      print "#ggu_hard_install $_";		# write old first
      print "use lib (\"$dir/fembin\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (perl): $_" if $::verbose;
      next;
    } else {
      print STDERR "$ARGV ($dir) (none): $_" if $::verbose;
    }
  }

  if( /^FEMDIR_INSTALL/ ) {
      $change = 1;
      print "SHYFEM_INSTALL=$dir        # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (shell_install): $_" if $::verbose;
  }

  print;
 }
}

sub reset_fem {

 while(<>) {

  if( /ggu_hard_install$/ ) {
    $change = 1;
    next;
  } elsif( s/^#ggu_hard_install // ) {
    $change = 1;
  }

  print;
 }
}

#-------------------------------------------------------------

