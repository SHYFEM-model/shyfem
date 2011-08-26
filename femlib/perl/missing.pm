#!/usr/bin/perl -w
#
# time series (TS) interpolation utilities
#
# example of usage:
#
# use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
#
# use missing;
#
# my $ms = new missing;
#
##############################################################
#
# $ms->set_flag($flag)			# example -999
# $ms->set_verbose($verbose)		# 1=verbose  0=not verbose
# $ms->interpolate($ts)			# ts is pointer to TS
# $ms->substitute($to,$from)		# substitute into $to from $from
#					# $to and $from are pointers to TS
#
##############################################################

use strict;

package missing;

##############################################################

sub new
{
    my $self;

    $self =     {
                         flag             =>      -900.
                        ,verbose          =>      1
                };

    bless $self;
    return $self;
}

###############################################################################

sub set_flag {

  my ($self,$flag) = @_;

  $self->{flag} = $flag;
}

sub set_verbose {

  my ($self,$verbose) = @_;

  $self->{verbose} = $verbose;
}

sub mark_as_bad {

  my ($self,$ts,$min,$max) = @_;

  my $n = @$ts;
  my $flag = $self->{flag};

  for(my $i=0;$i<$n;$i++) {
    my $val = $ts->[$i];
    if( $val < $min or $val > $max ) {
      $ts->[$i] = $flag;
      #print STDERR "MIN_MAX: $min $max $val\n";
    }
  }
}

sub interpolate {

  my ($self,$ts) = @_;

  my $n = @$ts;

  my $ivalid = -1;
  my $inblock = 0;
  my @flag = ();	#array to indicate where invalid data has been found
			#(may be ignored by calling program)

  for(my $i=0;$i<$n;$i++) {
    if( $self->is_valid($$ts[$i]) ) {	#valid
      if( $inblock ) {		#interpolate (end of invalid block found)
        my $istart = $ivalid;
        my $iend = $i;
        $self->intp_block($ts,$istart,$iend);
        $inblock = 0;
      }
      $ivalid = $i;
    } else {
     $inblock++;
     $flag[$i] = 1;
    }
  }

  if( $inblock ) {		#we are still in an invalid block -> finish
    if( $ivalid == -1 ) {
      die "Cannot interpolate: no good value\n";
    }
    $self->intp_block($ts,$ivalid,$n);
  }

  return \@flag;
}
  
sub intp_block {		#internal routine -> do not call directly

  my ($self,$ts,$i0,$i1) = @_;

  my $n = @$ts;
  my $verbose = $self->{verbose};

  my ($v0,$v1);

  if( $i0 < 0 ) {
    $v0 = $$ts[$i1];
  } else {
    $v0 = $$ts[$i0];
  }
  if( $i1 >= $n ) {
    $v1 = $$ts[$i0];
  } else {
    $v1 = $$ts[$i1];
  }
  
  my $aux = ($v1-$v0) / ($i1-$i0);

  for( my $i=$i0+1 ; $i<$i1 ; $i++ ) {
    my $val = $$ts[$i];
    my $newval = $v0 + ($i-$i0) * $aux;
    print STDERR "interpolating $i  $val -> $newval\n" if $verbose;
    $$ts[$i] = $newval;
  }
}  

sub substitute {

  my ($self,$to,$from) = @_;

  my $n = @$to - 1;
  my $verbose = $self->{verbose};
  my @flag = ();	#array to indicate where invalid data has been found
			#(may be ignored by calling program)

  for my $i (0..$n) {
    my $val = $$to[$i];
    unless( $self->is_valid($val) ) {
      my $newval = $$from[$i];
      print STDERR "substituting: $i  $val -> $newval\n" if $verbose;
      $$to[$i] = $newval;;
      $flag[$i] = 1;
    }
  }

  return \@flag;
}

sub is_valid {
 
    my ($self,$value) = @_;

    my $flag = $self->{flag};

    if( $value =~ /[a-z]/i ) {
      return 0;
    } elsif( $value > $flag ) {
      return 1;
    } else {
      return 0;
    }
}

################################
1;
################################

