c shypart.f, v 1.10 2016/02/16 CMCC 

c----------------------------------------------------------------

	program shypart_main

	use mod_geom    !ok
	use levels      !ok
	use basin       !ok
        use shypart     !ok

c----------------------------------------------------------------

c include files

	include 'param.h'
	include 'mkonst.h'
	include 'pkonst.h'
	include 'femtime.h'

c local variables

	logical bdebout
	integer iwhat,levdbg
	integer date,time
	integer nthreads
	integer*8 count1,count2, count_rate, count_max
	real time1,time2
	double precision timer

	real getpar

        integer ierr

	bdebout = .false.


	call cpu_time(time1)
	call system_clock(count1, count_rate, count_max)
!$      timer = omp_get_wtime() 

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%% code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call shyfem_copyright('3D FEM model')

c-----------------------------------------------------------
c dimensions
c-----------------------------------------------------------

c-----------------------------------------------------------
c read STR file
c-----------------------------------------------------------

	call cstinit

	call cstfile				!read STR and basin

        call shypart_init(.true.)

        call shypart_setup

        call partPHG

c        call deallocate_array

        call shypart_finalize

        end

c*****************************************************************

	subroutine debug_output(it)

	use mod_meteo
	use mod_waves
	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_gotm_aux
	use mod_ts
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin

	implicit none

	integer it

	!include 'param.h'

	write(66) it

	call debug_output_record(3*nel,3,hm3v,'hm3v')
	call debug_output_record(nkn,1,xgv,'xgv')
	call debug_output_record(nkn,1,ygv,'ygv')

	call debug_output_record(3*nel,3,zeov,'zeov')
	call debug_output_record(3*nel,3,zenv,'zenv')
	call debug_output_record(nkn,1,zov,'zov')
	call debug_output_record(nkn,1,znv,'znv')

	call debug_output_record(nlvdi*nel,nlvdi,utlov,'utlov')
	call debug_output_record(nlvdi*nel,nlvdi,vtlov,'vtlov')
	call debug_output_record(nlvdi*nel,nlvdi,utlnv,'utlnv')
	call debug_output_record(nlvdi*nel,nlvdi,vtlnv,'vtlnv')

        call debug_output_record(nlvdi*nkn,nlvdi,saltv,'saltv')
        call debug_output_record(nlvdi*nkn,nlvdi,tempv,'tempv')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,visv,'visv')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,wlov,'wlov')
	call debug_output_record((nlvdi+1)*nkn,nlvdi+1,wlnv,'wlnv')

	call debug_output_record(nkn,1,z0bk,'z0bk')
	call debug_output_record(nkn,1,tauxnv,'tauxnv')
	call debug_output_record(nkn,1,tauynv,'tauynv')

	call debug_output_record(nlvdi*nel,nlvdi,hdeov,'hdeov')
        call debug_output_record(nlvdi*nel,nlvdi,hdenv,'hdenv')
        call debug_output_record(nlvdi*nkn,nlvdi,hdkov,'hdkov')
        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'hdknv')

        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'shearf2')
        call debug_output_record(nlvdi*nkn,nlvdi,hdknv,'buoyf2')

        call debug_output_record(nlvdi*nel,nlvdi,fxv,'fxv')
        call debug_output_record(nlvdi*nel,nlvdi,fyv,'fyv')
        call debug_output_record(nlvdi*nel,nlvdi,wavefx,'wavefx')
        call debug_output_record(nlvdi*nel,nlvdi,wavefy,'wavefy')
        call debug_output_record(nel,1,rfricv,'rfricv')

        call debug_output_record(nlvdi*nkn,nlvdi,momentxv,'momentxv')
        call debug_output_record(nlvdi*nkn,nlvdi,momentyv,'momentyv')
        !call debug_output_record((nlvdi+1)*nkn,nlvdi+1,vts,'vts')

	write(66) 0,0

	end

c*****************************************************************

	subroutine debug_output_record(ntot,nfirst,val,text)
	implicit none
	integer ntot,nfirst
	real val(ntot)
	character*(*) text
	character*80 text1
	text1=text
	write(66) ntot,nfirst
	write(66) text1
	write(66) val
	end

c*****************************************************************

