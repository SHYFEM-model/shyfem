#!/usr/bin/perl
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# analyses common blocks (per file)

@modules = ( 
		 basin, topo, drywet, geom, hydro, femtim, constant
		,simul, aux, level, depth, tsc, bound
		,wind, friction, grd, matrix, float
		,mercator, semiimpl, flux, volume, param
		,nos, extra, qflx, namelist, graph, list
	   );

@basin = ( nkonst, descrr, nen3v, ipv, ipev, iarv, xgv, ygv, hm3v );
@topo = ( ilinkv, lenkv, kantv, ieltv, inodv, winv, dxv, dyv );
@drywet = ( iwegv );
@geom = ( ev );
@hydro = ( uov, unv, vov, vnv, zov, znv, zeov, zenv, up0v, vp0v, xv );
@femtim = ( femtim );
@constant = ( mkonst, pkonst );
@simul = ( descrp );
@aux = ( v1v, v2v, v3v, v4v, v5v, uedif, vedif );
#@dimension = ( dimdim );
@level = ( level, ilhv, ilhkv, hldv, hlv, hlhv );
@depth = ( hev, hkv );
@tsc = ( tnv, toev, snv, soev, cnv, coev, rtv, rsv, rcv, tempn, saltn, conzn );
@bound = ( bnd, irv, ierv, rzv, rqv, rlv, rhv, ruv, rvv, rrv, boundn );
@wind = ( tauxnv, tauynv, wxov, wxnv, wyov, wynv, ppv, pov, pnv 
		,winpar, wintim, windat );
@friction = ( czv, chezy, ausv );
@grd = ( vscale );
@matrix = ( amat, bmat );
@float = ( ffloat, tflfd, ieflv, iefdv, xpflv, ypflv, xpfdv, ypfdv );
@mercator = ( mercom, mcoo );
@semiimpl = ( semimi, semimr );
@flux = ( iflux, kfluxc );
@volume = ( ivol, kvolc );
@param = ( parsco, fnmsco );
@extra = ( extcom, knausc );
@nos = ( noscom, nosvar );
@qflx = ( qflxrn, qflxra, qflxro, qflxi );
@namelist = ( nrdcom );
@graph = ( ppp20 );
@list = ( lstloc );	#???

# read in common blocks (per file) -------------------------------

print "common blocks read...\n";

while(<>) {

  @f = split;

  $common = $f[0];
  $n = $f[1];
  $n =~ s/[\(\)]//g;
  $commonblock{$common} = $n;

  print "$n $common\n";
}

print "statistics of common blocks...\n";
&print_common_stats;

# prepare modules

print "available modules with common blocks...\n";

foreach $module (@modules) {
  $line = '@f = @' . "$module;";
  eval($line);	#now in @f are common blocks of module
  print "$module: ";
  &print_array(@f);
  foreach $item (@f) {
    $aux = "\/$item\/";
    unless( $commonblock{$aux} ) {
	die "error: common block not available: $aux\n";
    }
    $commonblock{$aux} = 0;
  }
}

print "statistics of common blocks...\n";
&print_common_stats;

##################################################

sub print_array {

  foreach $xxx (@_) {
    print "$xxx ";
  }
  print "\n";
}

##################################################

sub print_common_stats {

  $commax = 8;
  $commax = 20;
  $commax = 25;
  $commax = 40;
  $commax = 100;
  @nums = ();
  @coms = ();

  foreach $common (keys %commonblock) {
    $n = $commonblock{$common};
    $nums[$n]++;
    $coms[$n] .= "$common ";
  }

  $ntot = @nums;
  for($i=1;$i<=$ntot;$i++) {
    $n = $nums[$i];
    next unless $n;
    print "$i   $n";
    print "   $coms[$i]" if $n <= $commax;
    print "\n";
  }

}
