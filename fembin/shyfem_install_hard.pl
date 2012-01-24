#!/usr/bin/perl
#
# installs shyfem through hard references
#
# called by shell script
#
#---------------------------------------------

$reset = 0;
$dir = shift;
if( $dir eq "-reset" ) {
  $reset = 1;
  $dir = shift;
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
      print STDERR "$ARGV ($dir) (shell): $_";
    } elsif( /femlib\/perl/ ) {
      $change = 1;
      print;					# write old first
      print "use lib (\"$dir/femlib/perl\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (femlib/perl): $_";
      next;
    } elsif( /ENV/ and /fembin/ ) {
      $change = 1;
      print;					# write old first
      print "use lib (\"$dir/fembin\");     # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (perl): $_";
      next;
    } else {
      print STDERR "$ARGV ($dir) (none): $_";
    }
  }

  if( /^FEMDIR_INSTALL/ ) {
      $change = 1;
      print "SHYFEM_INSTALL=$dir        # ggu_hard_install\n";
      print STDERR "$ARGV ($dir) (shell_install): $_";
  }

  print;
 }
}

sub reset_fem {

 while(<>) {

  if( /ggu_hard_install$/ ) {
    $change = 1;
    next;
  }

  print;
 }
}

#-------------------------------------------------------------

