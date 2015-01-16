
!chezy.h		nardim			*
!bsig.h			nsidim, nkn		*
!bnd.h			nbndim,nbcdim		*G
!bfm_common.h		?			-
!bclfix.h		nkn,nel,nlv		C
!aquabc_aout.h		?			-
!aquabc.h		?			-
!param.h		-			-
!nlevel.h		-			-
!ts.h			nkn,nlv			F
!levels.h		nkn,nel,nlv		M
!meteo_aux.h		nkn			C
!meteo.h		nkn			C
!hydro.h		nkn,nel,nlv		M
!hydro_vel.h		nkn,nel,nlv		M
!hydro_print.h		nkn,nlv			M
!geom_dynamic.h		nkn,nel			H
!geom.h			nkn,nel,nlkdim		*G
!geom_aux.h		nkn,nel			G
!tides.h		nkn			C
!diff_visc_fric.h	nkn,nel,nlv		M
!waves.h		nkn,nlv			C
!bound_geom.h		nkn,nrbdim	-> subbndo.h		*G
!hydro_baro.h		nel			M
!area.h			nkn,nlv			F
!sinking.h		nkn,nlv			C
!simul.h		-			-
!aux_array.h		nkn,nel,nlv		F
!turbulence.h		nkn,nlv			C
!bound_dynamic.h	nkn,nlv			F
!bound_names.h		nbcdim		-> bnd.h, subboxa.h	*G
!diff_aux.h		nel			F
!bnd_aux.h		nkn,nel			H
!volcomp.h		nfxdim			*G
!nudging.h		nkn			H
!gotm_aux.h		nkn,nel,nlv		F
!conz.h			nkn,nlv,ncsdim		*G
!nohyd.h		nkn,nlv			H
!extra.h		nexdim			*
!depth.h		nkn,nel,nlv		H
!basin.h		nkn,nel			G
!roughness.h		nkn			M
!const_aux.h		nkn,nlv			F
!debug_aux1.h		-			-
!debug_aux2.h		-			-
!coords.h		-			-
!coords_gb.h		-			-
!coords_utm.h		-			-
!coords_cpp.h		-			-
!internal.h		nkn,nel,nlv		H
!sigma.h		-			-
!histo.h		ndim			-
!stab.h			ndim			-
!subgrd.h		-			-
!semi.h			-			-
!subnls.h		-			-
!fluidmud.h		nkn,nlv			C
!reg.h			-			-

!nfxdim:	flxinf.f splitflx.f subflxa.f volinf.f

! C 7    M 5    F 5    H 5    Z 5    G 8

! Chris: bclfix.h meteo_aux.h meteo.h tides.h waves.h sinking.h fluidmud.h
!		turbulence.h
! Marco: roughness.h levels.h hydro.h hydro_vel.h hydro_print.h hydro_baro.h
!		diff_visc_fric.h
! Fra: ts.h area.h aux_array.h bound_dynamic.h diff_aux.h
!		gotm_aux.h const_aux.h
! Michol: geom_dynamic.h bnd_aux.h nudging.h nohyd.h depth.h internal.h
!	subroutine fem_alloc
! Georg: bnd.h geom.h geom_aux.h (link.h) bound_geom.h volcomp.h
!	conz.h basin.h 

	subroutine text_alloc
	end
