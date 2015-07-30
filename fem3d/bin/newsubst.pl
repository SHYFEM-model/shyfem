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
$::write = 0 unless $::write;
$::debug = 0 unless $::debug;

$::subst = 0 unless $::subst;
$::clean = 0 unless $::clean;
$::include = 0 unless $::include;
$::inc2use = 0 unless $::inc2use;
$::revert = 0 unless $::revert;
$::check = 0 unless $::check;

#----------------------------------------------------------------

# newini.f		375: hlv(2)
#
# in sublnka.f subflx3d : links.h does not have nkndim etc...
# clean boundary names (???)
# subgotm: 214 -> comment shearf2 buoyf2 once gotm_aux.h is included

#----------------------------------------------------------------

$fortran->{debug} = $::debug;

$fortran->set_ignore_files(\@ignore_files);
$fortran->set_ignore_routines(\@ignore_routines);

$fortran->parse_files($::quiet);

#----------------------------------------------------------------

if( $::subst ) {
  subst_femtim($fortran);
  subst_konst($fortran);
  subst_ts($fortran);
  subst_levels($fortran);
  subst_depth($fortran);
  subst_meteo($fortran);
  subst_hydro($fortran);
  subst_hydro_vel($fortran);
  subst_hydro_print($fortran);
  subst_hydro_plot($fortran);
  subst_plot_aux($fortran);
  subst_plot_supout($fortran);
  subst_basin($fortran);
  subst_geom_dynamic($fortran);
  subst_geom($fortran);
  subst_tides($fortran);
  subst_diff_visc_fric($fortran);
  subst_waves($fortran);
  subst_bound_geom($fortran);
  subst_hydro_baro($fortran);
  subst_sinking($fortran);
  subst_simul($fortran);
  subst_area($fortran);
  subst_aux_array($fortran);
  subst_turbulence($fortran);
  subst_bound_dynamic($fortran);
  subst_bound_names($fortran);
  subst_roughness($fortran);
  subst_internal($fortran);
  subst_diff_aux($fortran);
  subst_bnd_aux($fortran);
  subst_fluidmud($fortran);
  subst_volcomp($fortran);
  subst_nudging($fortran);
  subst_gotm($fortran);
  subst_conz($fortran);
  subst_nohyd($fortran);
  subst_extra($fortran);
  subst_const_aux($fortran);
  subst_debug($fortran);
  subst_coords($fortran);
  subst_sigma($fortran);
  subst_histo($fortran);
  subst_stab($fortran);
  subst_grd($fortran);
  subst_semi($fortran);
  subst_subnls($fortran);
  subst_reg($fortran);

  treat_includes($fortran);
  #check_common($fortran);
} elsif( $::clean ) {
  clean_files($fortran);
} elsif( $::include ) {
  treat_includes($fortran);
} elsif( $::inc2use ) {
  inc2use($fortran);
} elsif( $::revert ) {
  treat_reverts($fortran);
} elsif( $::check ) {
  check_common($fortran);
}

$fortran->write_files($::quiet) if $::write;

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub subst_reg {
  my $fortran = shift;
  subst_common($fortran,"ppp20","reg.h");
}

sub subst_subnls {
  my $fortran = shift;
  subst_common($fortran,"secsec","subnls.h");
  subst_common($fortran,"nrdcom","subnls.h");
}

sub subst_semi {
  my $fortran = shift;
  subst_common($fortran,"semimi","semi.h");
  subst_common($fortran,"semimr","semi.h");
}

sub subst_grd {
  my $fortran = shift;
  subst_common($fortran,"vscale","subgrd.h");
  subst_common($fortran,"grdcom_i","subgrd.h");
  subst_common($fortran,"grdcom_r","subgrd.h");
  subst_common($fortran,"grdcom_c","subgrd.h");
}

sub subst_stab {
  my $fortran = shift;
  subst_common($fortran,"istab","stab.h");
  subst_common($fortran,"rstab","stab.h");
}

sub subst_histo {
  my $fortran = shift;
  subst_common($fortran,"ihisto","histo.h");
  subst_common($fortran,"rhisto","histo.h");
}

sub subst_sigma {
  my $fortran = shift;
  subst_common($fortran,"nsigma_com","sigma.h");
  subst_common($fortran,"hsigma_com","sigma.h");
}

sub subst_coords {
  my $fortran = shift;
  subst_common($fortran,"coords1","coords.h");
  subst_common($fortran,"proj_param01","coords.h");
  subst_common($fortran,"proj_param02","coords.h");
  subst_common($fortran,"proj_param11","coords.h");
  subst_common($fortran,"proj_param12","coords.h");
  subst_common($fortran,"proj_param13","coords.h");
  subst_common($fortran,"proj_param21","coords.h");
  subst_common($fortran,"gb_param1","coords_gb.h");
  subst_common($fortran,"gb_param2","coords_gb.h");
  subst_common($fortran,"utm_param1","coords_utm.h");
  subst_common($fortran,"utm_param2","coords_utm.h");
  subst_common($fortran,"utm_param3","coords_utm.h");
  subst_common($fortran,"cpp_param1","coords_cpp.h");
}

sub subst_debug {
  my $fortran = shift;
  subst_common($fortran,"comdebug","debug_aux1.h");
  subst_common($fortran,"debug_ggu","debug_aux2.h");
}

sub subst_extra {
  my $fortran = shift;
  subst_common($fortran,"knausc","param_dummy.h","extra.h");
}

sub subst_nohyd {
  my $fortran = shift;
  subst_common($fortran,"qpnv","param_dummy.h","nohyd.h");
  subst_common($fortran,"qpov","param_dummy.h","nohyd.h");
}

sub subst_conz {
  my $fortran = shift;
  subst_common($fortran,"conzv","param_dummy.h","conz.h");
  subst_common($fortran,"cnv","param_dummy.h","conz.h");
}

sub subst_gotm {
  my $fortran = shift;
  subst_common($fortran,"numv_gotm","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"nuhv_gotm","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"tken_gotm","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"eps_gotm","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"rls_gotm","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"shearf2","param_dummy.h","gotm_aux.h");
  subst_common($fortran,"buoyf2","param_dummy.h","gotm_aux.h");
}

sub subst_nudging {
  my $fortran = shift;
  subst_common($fortran,"andgzv","param_dummy.h","nudging.h");
}

sub subst_volcomp {
  my $fortran = shift;
  subst_common($fortran,"kvolc","param_dummy.h","volcomp.h");
  subst_common($fortran,"ivol","param_dummy.h","volcomp.h");
}

sub subst_fluidmud {
  my $fortran = shift;
  subst_common($fortran,"rhosed","param_dummy.h","fluidmud.h");
  subst_common($fortran,"dm0","param_dummy.h","fluidmud.h");
  subst_common($fortran,"nf","param_dummy.h","fluidmud.h");
  subst_common($fortran,"z0bkmud","param_dummy.h","fluidmud.h");
  subst_common($fortran,"mudc","param_dummy.h","fluidmud.h");
  subst_common($fortran,"rhomud","param_dummy.h","fluidmud.h");
  subst_common($fortran,"visv_yield","param_dummy.h","fluidmud.h");
  subst_common($fortran,"difv_yield","param_dummy.h","fluidmud.h");
  subst_common($fortran,"lambda","param_dummy.h","fluidmud.h");
  subst_common($fortran,"vts","param_dummy.h","fluidmud.h");
  subst_common($fortran,"dmf_mud","param_dummy.h","fluidmud.h");
  subst_common($fortran,"wprvs","param_dummy.h","fluidmud.h");
}

sub subst_bnd_aux {
  my $fortran = shift;
  subst_common($fortran,"ruv","param_dummy.h","bnd_aux.h");
  subst_common($fortran,"rvv","param_dummy.h","bnd_aux.h");
  subst_common($fortran,"crad","param_dummy.h","bnd_aux.h");
}

sub subst_internal {
  my $fortran = shift;
  subst_common($fortran,"rdistv","param_dummy.h","internal.h");
  subst_common($fortran,"fcorv","param_dummy.h","internal.h");
  subst_common($fortran,"fxv","param_dummy.h","internal.h");
  subst_common($fortran,"fyv","param_dummy.h","internal.h");
  subst_common($fortran,"iuvfix","param_dummy.h","internal.h");
  subst_common($fortran,"ddxv","param_dummy.h","internal.h");
  subst_common($fortran,"ddyv","param_dummy.h","internal.h");
}

sub subst_diff_aux {
  my $fortran = shift;
  subst_common($fortran,"wdifhv","param_dummy.h","diff_aux.h");
}

sub subst_const_aux {
  my $fortran = shift;
  subst_common($fortran,"const3d","param_dummy.h","const_aux.h");
}

sub subst_bound_names {
  my $fortran = shift;
  subst_common($fortran,"boundn","param_dummy.h","bound_names.h");
  subst_common($fortran,"conzn","param_dummy.h","bound_names.h");
  subst_common($fortran,"saltn","param_dummy.h","bound_names.h");
  subst_common($fortran,"tempn","param_dummy.h","bound_names.h");
  subst_common($fortran,"bio2dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"sed2dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"mud2dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"lam2dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"dmf2dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"tox3dn","param_dummy.h","bound_names.h");
  subst_common($fortran,"bfm1bc","param_dummy.h","bound_names.h");
  subst_common($fortran,"bfm2bc","param_dummy.h","bound_names.h");
  subst_common($fortran,"bfm3bc","param_dummy.h","bound_names.h");
  subst_common($fortran,"vel3dn","param_dummy.h","bound_names.h");
}

sub subst_bound_dynamic {
  my $fortran = shift;
  subst_common($fortran,"mfluxv","param_dummy.h","bound_dynamic.h");
  subst_common($fortran,"rqpsv","param_dummy.h","bound_dynamic.h");
  subst_common($fortran,"rqdsv","param_dummy.h","bound_dynamic.h");
  subst_common($fortran,"rzv","param_dummy.h","bound_dynamic.h");
  subst_common($fortran,"rqv","param_dummy.h","bound_dynamic.h");
}

sub subst_turbulence {
  my $fortran = shift;
  subst_common($fortran,"tken","param_dummy.h","turbulence.h");
  subst_common($fortran,"eps","param_dummy.h","turbulence.h");
  subst_common($fortran,"rls","param_dummy.h","turbulence.h");
}

sub subst_roughness {
  my $fortran = shift;
  subst_common($fortran,"z0bk","param_dummy.h","roughness.h");
  subst_common($fortran,"z0sk","param_dummy.h","roughness.h");
}

sub subst_aux_array {
  my $fortran = shift;
  subst_common($fortran,"v1v","param_dummy.h","aux_array.h");
  subst_common($fortran,"v2v","param_dummy.h","aux_array.h");
  subst_common($fortran,"v3v","param_dummy.h","aux_array.h");
  subst_common($fortran,"ve1v","param_dummy.h","aux_array.h");
  subst_common($fortran,"saux1","param_dummy.h","aux_array.h");
  subst_common($fortran,"saux2","param_dummy.h","aux_array.h");
  subst_common($fortran,"saux3","param_dummy.h","aux_array.h");
  subst_common($fortran,"saux4","param_dummy.h","aux_array.h");
  subst_common($fortran,"sauxe1","param_dummy.h","aux_array.h");
  subst_common($fortran,"sauxe2","param_dummy.h","aux_array.h");
}

sub subst_area {
  my $fortran = shift;
  subst_common($fortran,"areakv","param_dummy.h","area.h");
}

sub subst_sinking {
  my $fortran = shift;
  subst_common($fortran,"wsinkv","param_dummy.h","sinking.h");
}

sub subst_simul {
  my $fortran = shift;
  subst_common($fortran,"descrp","param_dummy.h","simul.h");
}

sub subst_hydro_baro {
  my $fortran = shift;
  subst_common($fortran,"uov","param_dummy.h","hydro_baro.h");
  subst_common($fortran,"vov","param_dummy.h","hydro_baro.h");
  subst_common($fortran,"unv","param_dummy.h","hydro_baro.h");
  subst_common($fortran,"vnv","param_dummy.h","hydro_baro.h");
}

sub subst_bound_geom {
  my $fortran = shift;
  subst_common($fortran,"ierv","param_dummy.h","bound_geom.h");
  subst_common($fortran,"irv","param_dummy.h","bound_geom.h");
  subst_common($fortran,"iopbnd","param_dummy.h","bound_geom.h");
  subst_common($fortran,"rhv","param_dummy.h","bound_geom.h");
  subst_common($fortran,"rlv","param_dummy.h","bound_geom.h");
  subst_common($fortran,"rrv","param_dummy.h","bound_geom.h");
}

sub subst_waves {
  my $fortran = shift;
  subst_common($fortran,"waveh","param_dummy.h","waves.h");
  subst_common($fortran,"wavep","param_dummy.h","waves.h");
  subst_common($fortran,"wavepp","param_dummy.h","waves.h");
  subst_common($fortran,"waved","param_dummy.h","waves.h");
  subst_common($fortran,"waveov","param_dummy.h","waves.h");
  subst_common($fortran,"wavefx","param_dummy.h","waves.h");
  subst_common($fortran,"wavefy","param_dummy.h","waves.h");
}

sub subst_diff_visc_fric {
  my $fortran = shift;
  subst_common($fortran,"rfricv","param_dummy.h","diff_visc_fric.h");
  subst_common($fortran,"czv","param_dummy.h","diff_visc_fric.h");
  subst_common($fortran,"austv","param_dummy.h","diff_visc_fric.h");
  subst_common($fortran,"difhv","param_dummy.h","diff_visc_fric.h");
  #subst_common($fortran,"wdifhv","param_dummy.h","diff_visc_fric.h");
  subst_common($fortran,"visv","param_dummy.h","diff_visc_fric.h");
  subst_common($fortran,"difv","param_dummy.h","diff_visc_fric.h");
}

sub subst_tides {
  my $fortran = shift;
  subst_common($fortran,"xgeov","param_dummy.h","tides.h");
  subst_common($fortran,"ygeov","param_dummy.h","tides.h");
  subst_common($fortran,"zeqv","param_dummy.h","tides.h");
}

sub subst_geom {
  my $fortran = shift;
  subst_common($fortran,"ilinkv","param_dummy.h","geom.h");
  subst_common($fortran,"lenkv","param_dummy.h","geom.h");
  subst_common($fortran,"linkv","param_dummy.h","geom.h");
  subst_common($fortran,"ieltv","param_dummy.h","geom.h");
  subst_common($fortran,"kantv","param_dummy.h","geom.h");
  subst_common($fortran,"dxv","param_dummy.h","geom.h");
  subst_common($fortran,"dyv","param_dummy.h","geom.h");
}

sub subst_geom_dynamic {
  my $fortran = shift;
  subst_common($fortran,"iwegv","param_dummy.h","geom_dynamic.h");
  subst_common($fortran,"iwetv","param_dummy.h","geom_dynamic.h");
  subst_common($fortran,"inodv","param_dummy.h","geom_dynamic.h");
}

sub subst_basin {
  my $fortran = shift;
  subst_common($fortran,"descrr","param_dummy.h","basin.h");
  subst_common($fortran,"xgv","param_dummy.h","basin.h");
  subst_common($fortran,"ygv","param_dummy.h","basin.h");
  subst_common($fortran,"hm3v","param_dummy.h","basin.h");
  subst_common($fortran,"nen3v","param_dummy.h","basin.h");
  subst_common($fortran,"ipev","param_dummy.h","basin.h");
  subst_common($fortran,"ipv","param_dummy.h","basin.h");
  subst_common($fortran,"iarv","param_dummy.h","basin.h");
  subst_common($fortran,"iarnv","param_dummy.h","basin.h");
}

sub subst_plot_supout {
  my $fortran = shift;
  subst_common($fortran,"eoseos","supout.h");
  subst_common($fortran,"nosnos","supout.h");
  subst_common($fortran,"wavwav","supout.h");
  subst_common($fortran,"femfem","supout.h");
  subst_common($fortran,"fvlfvl","supout.h");
  subst_common($fortran,"ousous","supout.h");
}

sub subst_plot_aux {
  my $fortran = shift;
  subst_common($fortran,"fvlv","param_dummy.h","plot_aux.h");
  subst_common($fortran,"arfvlv","param_dummy.h","plot_aux.h");
  subst_common($fortran,"wauxv","param_dummy.h","plot_aux.h");
  subst_common($fortran,"bwater","param_dummy.h","plot_aux.h");
  subst_common($fortran,"bkwater","param_dummy.h","plot_aux.h");
  subst_common($fortran,"hetv","param_dummy.h","plot_aux.h");
  subst_common($fortran,"het3v","param_dummy.h","plot_aux.h");
  #subst_common($fortran,"hl","param_dummy.h","plot_aux.h");
  subst_common($fortran,"parray","param_dummy.h","plot_aux.h");
  subst_common($fortran,"p3","param_dummy.h","plot_aux.h");
}

sub subst_hydro_plot {
  my $fortran = shift;
  subst_common($fortran,"uv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"vv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"uvnv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"vvnv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"usnv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"vsnv","param_dummy.h","hydro_plot.h");
  subst_common($fortran,"wsnv","param_dummy.h","hydro_plot.h");
}

sub subst_hydro_print {
  my $fortran = shift;
  subst_common($fortran,"uprv","param_dummy.h","hydro_print.h");
  subst_common($fortran,"vprv","param_dummy.h","hydro_print.h");
  subst_common($fortran,"upro","param_dummy.h","hydro_print.h");
  subst_common($fortran,"vpro","param_dummy.h","hydro_print.h");
  subst_common($fortran,"wprv","param_dummy.h","hydro_print.h");
  subst_common($fortran,"up0v","param_dummy.h","hydro_print.h");
  subst_common($fortran,"vp0v","param_dummy.h","hydro_print.h");
  subst_common($fortran,"xv","param_dummy.h","hydro_print.h");
}

sub subst_hydro_vel {
  my $fortran = shift;
  subst_common($fortran,"ulov","param_dummy.h","hydro_vel.h");
  subst_common($fortran,"ulnv","param_dummy.h","hydro_vel.h");
  subst_common($fortran,"vlov","param_dummy.h","hydro_vel.h");
  subst_common($fortran,"vlnv","param_dummy.h","hydro_vel.h");
  subst_common($fortran,"wlov","param_dummy.h","hydro_vel.h");
  subst_common($fortran,"wlnv","param_dummy.h","hydro_vel.h");
}

sub subst_hydro {
  my $fortran = shift;
  subst_common($fortran,"zov","param_dummy.h","hydro.h");
  subst_common($fortran,"znv","param_dummy.h","hydro.h");
  subst_common($fortran,"zeov","param_dummy.h","hydro.h");
  subst_common($fortran,"zenv","param_dummy.h","hydro.h");
  subst_common($fortran,"utlov","param_dummy.h","hydro.h");
  subst_common($fortran,"utlnv","param_dummy.h","hydro.h");
  subst_common($fortran,"vtlov","param_dummy.h","hydro.h");
  subst_common($fortran,"vtlnv","param_dummy.h","hydro.h");
}

sub subst_meteo {
  my $fortran = shift;
  subst_common($fortran,"metwbt","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"metws","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"metrain","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"ppv","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"metrad","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"methum","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"mettair","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"metcc","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"tauxnv","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"tauynv","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"wxv","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"wyv","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"windcd","param_dummy.h","meteo_aux.h");
  subst_common($fortran,"evapv","param_dummy.h","meteo_aux.h");
}

sub subst_depth {
  my $fortran = shift;
  subst_common($fortran,"hkv","param_dummy.h","depth.h");
  subst_common($fortran,"hev","param_dummy.h","depth.h");
  subst_common($fortran,"hdknv","param_dummy.h","depth.h");
  subst_common($fortran,"hdkov","param_dummy.h","depth.h");
  subst_common($fortran,"hdenv","param_dummy.h","depth.h");
  subst_common($fortran,"hdeov","param_dummy.h","depth.h");
  subst_common($fortran,"hkv_min","param_dummy.h","depth.h");
  subst_common($fortran,"hkv_max","param_dummy.h","depth.h");
}

sub subst_levels {
  my $fortran = shift;
  subst_common($fortran,"ilhkv","param_dummy.h","levels.h");
  subst_common($fortran,"ilhv","param_dummy.h","levels.h");
  subst_common($fortran,"hlv","param_dummy.h","levels.h");
  subst_common($fortran,"hldv","param_dummy.h","levels.h");
  subst_common($fortran,"ilmv","param_dummy.h","levels.h");
  subst_common($fortran,"ilmkv","param_dummy.h","levels.h");
}

sub subst_ts {
  my $fortran = shift;
  subst_common($fortran,"rhov","ts.h","param_dummy.h");
  subst_common($fortran,"saltv","ts.h","param_dummy.h");
  subst_common($fortran,"tempv","ts.h","param_dummy.h");
  subst_common($fortran,"sobsv","ts.h","param_dummy.h");
  subst_common($fortran,"tobsv","ts.h","param_dummy.h");
  subst_common($fortran,"rtauv","ts.h","param_dummy.h");
  subst_common($fortran,"bpresv","ts.h","param_dummy.h");
  subst_common($fortran,"bpresxv","ts.h","param_dummy.h");
  subst_common($fortran,"bpresyv","ts.h","param_dummy.h");
}

sub subst_konst {
  my $fortran = shift;
  subst_common($fortran,"mkonst","mkonst.h");
  subst_common($fortran,"pkonst","pkonst.h");
  subst_common($fortran,"nkonst","nbasin.h");
  subst_common($fortran,"level","nlevel.h");
}
sub subst_femtim {
  my $fortran = shift;
  subst_common($fortran,"femtim","femtime.h");
  subst_common($fortran,"femtimu","femtime.h");
}

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

sub check_common {

  my ($fortran) = @_;

  my %count = ();
  my %files = ();

  my @list = sort keys %{$fortran->{all_routines}};
  foreach my $rname (@list) {
    my $ritem = $fortran->{all_routines}->{$rname};
    my @clist = keys %{$ritem->{common}};
    my $file = $ritem->{file};
    foreach my $common (@clist) {
      $count{$common}++;
      $files{$common} .= "$file ";
    }
  }

  #my @sorted = sort keys %count;
  my @sorted = sort { $count{$a} <=> $count{$b} } keys %count;

  my $total = 0;
  my $number = 0;
  foreach my $common (@sorted) {
    my $files = condense_list($files{$common});
    my $nfiles = @$files;
    my $count = $count{$common};
    $total += $count;
    $number++;
    print "$count    ($nfiles)   $common\n";
  }
  print "total: $total  $number\n";
  #return;

  print "more than one file:\n";
  foreach my $common (@sorted) {
    my $files = condense_list($files{$common});
    my $nfiles = @$files;
    next if $nfiles <= 1;
    my $count = $count{$common};
    print "$count    ($nfiles)   $common\n";
  }

  my %byfile = ();
  print "one file:\n";
  foreach my $common (@sorted) {
    my $files = condense_list($files{$common});
    my $nfiles = @$files;
    next if $nfiles > 1;
    my $count = $count{$common};
    my $file = $files->[0];
    $byfile{$file} .= "$common ";
    print "$count    ($file)   $common\n";
  }

  print "by file:\n";
  foreach my $file (keys %byfile) {
    my $common = $byfile{$file};
    my $commons = condense_list($common);
    my $ncommons = @$commons;
    print "$file ($ncommons):    $common\n";
  }
}

sub condense_list {

  my $files = shift;

  my @files = split(/\s+/,$files);
  my %aux = ();
  foreach (@files) {
    $aux{$_}++;
  }
  my @aux = sort keys %aux;
  return \@aux
}

sub copy_include {

  my @include = @_;

  my $hdir = "/home/georg/fem/fem3d/common_h";

  foreach my $include (@include) {
    my $file = "$hdir/$include";
    if( -f $include ) {
      print STDERR "include file $include existing... not copying\n";
    } elsif( not -f $file ) {
      die "no include file $include available... aborting\n";
    } else {
      print STDERR "copy include file $include\n";
      `cp $file .`
    }
  }
}

sub subst_common {

  my ($fortran,$common,@include) = @_;

  my $include = join(" ",@include);
  print STDERR "  substituting common /$common/ with include $include\n";

  my @list = sort keys %{$fortran->{all_routines}};
  foreach my $rname (@list) {
    my $ritem = $fortran->{all_routines}->{$rname};
    my $name = $ritem->{name};
    my $file = $ritem->{file};

    if( $ritem->{common}->{$common} ) {
      print STDERR "    routine $name contains common... substituting\n";
      substitute_common($fortran,$ritem,$common,@include);
      copy_include(@include);		# copy only if changed
      $fortran->set_changed($file);
    }
  }
}

sub substitute_common {

  my ($fortran,$ritem,$common,@include) = @_;

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
	      handle_include($fortran,$ritem,\@new,\@include);
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

sub handle_include {

  my ($fortran,$ritem,$new,$include) = @_;

  #my $line = join(" ",@$include);
  #print "treat include: $line\n";

  foreach my $inc (@$include) {
    if( is_include_compatible($ritem,$inc) ) {
      push(@$new,"\tinclude '$inc' !COMMON_GGU_SUBST");
      $ritem->{include}->{$inc}++;
    }
  }
}

sub is_include_compatible {

  my ($ritem,$inc) = @_;

  return 0 if $ritem->{include}->{$inc};

  if( $inc eq "param_dummy.h" ) {
    if( $ritem->{include}->{"param.h"} ) {
      return 0;
    }
  }

  return 1;
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
  $line =~ s/^\*\(?\d+\)?// if $type eq "character";

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
#--------------------------------------------------------------
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

  write_newfile($filename,".clean",\@new);
}

sub write_newfile {

  my ($filename,$post,$new) = @_;

  my $newfile = $filename . $post;
  print STDERR "writing file $filename to $newfile\n";
  open(NEW,">$newfile");
  foreach my $line (@$new) {
    print NEW "$line\n";
  }
  close(NEW);
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub treat_reverts {

  my ($fortran) = @_;

  foreach my $filename (keys %{$fortran->{files}}) {
    treat_revert($fortran,$filename)
  }
}

sub treat_revert {

  my ($fortran,$filename) = @_;

  my $fitem = $fortran->{files}->{$filename};
  my $sequence = $fitem->{sequence};

  my $changed = 0;
  my @new = ();
  foreach my $ritem (@$sequence) {
    my $code = $ritem->{code};
    my $name = $ritem->{name};
    foreach my $line (@$code) {
      if( $line =~ s/^COMMON_GGU_DELETED// ) { $changed++; }
      if( $line =~ /\s*!COMMON_GGU_SUBST\s*$/ ) { $changed++; next; }
      push(@new,$line);
    }
  }

  return unless $changed;

  write_newfile($filename,".new",\@new);
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub treat_includes {

  my ($fortran) = @_;

  foreach my $filename (keys %{$fortran->{files}}) {
    treat_include($fortran,$filename)
  }
}

# tutto al primo include, subito dopo implicit none or dopo subroutine
# dopo use
# eliminare param - dummy
# cambiare dummy to param if param in file

sub treat_include {

  my ($fortran,$filename) = @_;

  my $fitem = $fortran->{files}->{$filename};
  my $sequence = $fitem->{sequence};

  my $param = 0;
  my $param_dummy = 0;
  foreach my $ritem (@$sequence) {
    $param++ if $ritem->{include}->{"param.h"};
    $param_dummy++ if $ritem->{include}->{"param_dummy.h"};
  }

  my $changed = 0;
  my @new = ();
  foreach my $ritem (@$sequence) {
    my $code = $ritem->{code};
    my $name = $ritem->{name};
    if( $ritem->{include}->{"param.h"} and
		$ritem->{include}->{"param_dummy.h"} ) {
      $changed += delete_include($ritem,"param_dummy.h");
    }
    if( $param and $ritem->{include}->{"param_dummy.h"} ) {
      $changed += substitute_include($ritem,"param_dummy.h","param.h");
    }
    if( $ritem->{include}->{"nbasin.h"} and
		$ritem->{include}->{"basin.h"} ) {
      $changed += delete_include($ritem,"nbasin.h");
    }
    if( $ritem->{include}->{"links.h"} and
		$ritem->{include}->{"geom.h"} ) {
      $changed += substitute_include($ritem,"geom.h","geom_aux.h");
    }

    my $line = "";
    my $ninc = -1;
    if( $ritem->{include}->{"param.h"} ) {
      ($ninc,$line) = get_include($ritem,"param.h");
    } elsif( $ritem->{include}->{"param_dummy.h"} ) {
      ($ninc,$line) = get_include($ritem,"param_dummy.h");
    }
    die "cannot find include (internal error)\n" if $ninc == 0;
    if( $ninc > 1 ) {
      $changed += make_first_include($ritem,$line);
    }

    $fortran->set_changed($filename) if $changed;
  }
}

sub get_include {

  my ($ritem,$inc) = @_;

  my $code = $ritem->{code};
  my $ninc = 0;

  foreach my $line (@$code) {
    my $include = is_include($line);
    if( $include ) {
      $ninc++;
      return ($ninc,$line) if( $include eq $inc );
    }
  }
  return (0,"");
}

sub substitute_include {

  my ($ritem,$inc,$subst) = @_;

  my $code = $ritem->{code};
  my @new = ();
  my $changed = 0;

  foreach my $line (@$code) {
    my $include = is_include($line);
    if( $include eq $inc ) {
      $line = "COMMON_GGU_DELETED$line";
      delete $ritem->{include}->{$inc};
      push(@new,$line);
      $line = "\tinclude '$subst' !COMMON_GGU_SUBST";
      $ritem->{include}->{$subst}++;
      $changed++;
    }
    push(@new,$line);
  }
  $ritem->{code} = \@new if $changed;
  return $changed;
}

sub delete_include {

  my ($ritem,$inc) = @_;

  my $code = $ritem->{code};
  my @new = ();
  my $changed = 0;

  foreach my $line (@$code) {
    my $include = is_include($line);
    if( $include eq $inc ) {
      $line = "COMMON_GGU_DELETED$line";
      delete $ritem->{include}->{$inc};
      $changed++;
    }
    push(@new,$line);
  }
  $ritem->{code} = \@new if $changed;
  return $changed;
}

sub make_first_include {

  my ($ritem,$pline) = @_;

  my $code = $ritem->{code};
  my @new = ();
  my $changed = 0;
  my $ninc = 0;

  foreach my $line (@$code) {
    if( $line eq $pline ) {
      $line = "COMMON_GGU_DELETED$pline";
      $changed++;
    }
    my $include = is_include($line);
    $ninc++ if $include;
    if( $include and $ninc == 1 ) {
      my $aline = "$pline !COMMON_GGU_SUBST";
      push(@new,$aline);
    }
    push(@new,$line);
  }
  $ritem->{code} = \@new if $changed;
  return $changed;
}

sub is_include {

  my $line = shift;

  $line =~ s/\s+//g;

  if( $line =~ /^include\'(\S+)\'/i or $line =~ /^include\"(\S+)\"/i ) {
    return $1;
  } else {
    return "";
  }
}

#--------------------------------------------------------------

sub include2use {

  my ($fortran,$include,$use,$only) = @_;

  my @list = sort keys %{$fortran->{all_routines}};
  foreach my $rname (@list) {
    my $ritem = $fortran->{all_routines}->{$rname};
    my $name = $ritem->{name};
    my $file = $ritem->{file};

    if( $ritem->{include}->{$include} ) {
      if( not $ritem->{use}->{$use} ) {		#not yet included
        print STDERR "    routine $name contains include... substituting\n";
        substitute_use($fortran,$ritem,$include,$use,$only);
        $ritem->{use}->{$use} = 1;
        $fortran->set_changed($file);
      } else {
	delete_include($ritem,$include);
        $fortran->set_changed($file);
      }
    }
  }
}

sub substitute_use {

  my ($fortran,$ritem,$include,$use,$only) = @_;

  my $code = $ritem->{code};
  my @new = ();
  my $nline = 0;
  my $used = 0;

  $only = ", only : $only" if $only;
  $only = "" unless $only;

  foreach my $l (@$code) {
    $nline++;
    if( $nline > 1 and not is_conti($l) ) {
      my $line = $fortran->clean_line($l);
      if( $line ) {				# first good line found
        $line =~ s/\s+//g;
	unless( $used ) {
          my $aline = "\tuse $use$only !COMMON_GGU_SUBST";
          push(@new,$aline);
	  unless( is_use($line) ) {		# add empty line
            push(@new,"");
          }
	  $used++;
	}
	if( is_include($line) eq $include ) {
	  $l = "COMMON_GGU_DELETED$l";
        }
      }
    }
    push(@new,$l);
  }

  $ritem->{code} = \@new;
}

sub is_use {

  my $line = shift;

  $line =~ s/\s+//g;

  if( $line =~ /^use(\w+)\b/i ) {
    return $1;
  } else {
    return "";
  }
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

sub inc2use {

  my $fortran = shift;

  include2use($fortran,"basin.h","basin");
  include2use($fortran,"nbasin.h","basin","nkn,nel,ngr,mbw");
  include2use($fortran,"levels.h","levels");
  include2use($fortran,"nlevel.h","levels","nlvdi,nlv");
  include2use($fortran,"ev.h","evgeom");
  include2use($fortran,"evmain.h","evgeom");

#  return;

  include2use($fortran,"hydro.h","mod_hydro");
  include2use($fortran,"hydro_vel.h","mod_hydro_vel");
  include2use($fortran,"hydro_print.h","mod_hydro_print");
  include2use($fortran,"hydro_baro.h","mod_hydro_baro");

  include2use($fortran,"diff_visc_fric.h","mod_diff_visc_fric");
  include2use($fortran,"roughness.h","mod_roughness");

#  return;

  include2use($fortran,"ts.h","mod_ts");
  include2use($fortran,"area.h","mod_area");
  include2use($fortran,"aux_array.h","mod_aux_array");
  include2use($fortran,"bound_dynamic.h","mod_bound_dynamic");
  include2use($fortran,"diff_aux.h","mod_diff_aux");
  include2use($fortran,"gotm_aux.h","mod_gotm_aux");

  include2use($fortran,"bnd_aux.h","mod_bnd_aux");
  include2use($fortran,"depth.h","mod_depth");
  include2use($fortran,"geom_dynamic.h","mod_geom_dynamic");
  include2use($fortran,"internal.h","mod_internal");
  include2use($fortran,"nohyd.h","mod_nohyd");
  include2use($fortran,"nudging.h","mod_nudging");

#  return;

  include2use($fortran,"bclfix.h","mod_bclfix");
  include2use($fortran,"fluidmud.h","mod_fluidmud");
  include2use($fortran,"sinking.h","mod_sinking");
  include2use($fortran,"turbulence.h","mod_turbulence");
  include2use($fortran,"waves.h","mod_waves");

  include2use($fortran,"meteo.h","mod_meteo");
  include2use($fortran,"meteo_aux.h","mod_meteo");

#  return;

  include2use($fortran,"geom.h","mod_geom");
  include2use($fortran,"geom_aux.h","mod_geom");
  include2use($fortran,"links.h","mod_geom");

  include2use($fortran,"conz.h","mod_conz");

  include2use($fortran,"bound_geom.h","mod_bound_geom");
  include2use($fortran,"testbndo.h","mod_bound_geom");
  include2use($fortran,"bnd.h","mod_bnd");
  include2use($fortran,"nbound.h","mod_bnd");
  include2use($fortran,"nbvdim.h","mod_bnd");

#  return;

  include2use($fortran,"tvd.h","mod_tvd");
  include2use($fortran,"tides.h","mod_tides");

  include2use($fortran,"subbndo.h","mod_bndo");
  include2use($fortran,"nudge.h","mod_nudge");

#  return;

  include2use($fortran,"plot_aux_3d.h","mod_plot3d");
  include2use($fortran,"plot_aux.h","mod_plot2d");
  include2use($fortran,"hydro_plot.h","mod_hydro_plot");

  include2use($fortran,"common.h","mod_system");
  include2use($fortran,"common_amat.h","mod_system");

  include2use($fortran,"grade.h","mod_adj_grade");
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

