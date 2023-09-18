!
! $Id: atoxi3d.f,v 1.6 2008-10-10 09:29:54 georg Exp $
!
! atoxi3d - Toxical routines from ARPAV - shell
!
! contents :
!
! revision log :
!
! 15.02.2006    ggu&fdp new routine atoxi3d for ARPAV (from bio3d)
! 23.03.2006    ggu     ntot eliminated
! 23.03.2006    ggu     changed time step to double precision
! 17.04.2008    ggu     new open boundary conditions 
! 22.04.2008    ggu     advection parallelized
! 23.04.2008    ggu     call to bnds_set_def() changed
! 09.10.2008    ggu     new call to confop
! 21.10.2014    ggu     converted to new boundary treatment
!
! notes :
!
! State variables used: (Wasp)
!
! todo :
!
! - * FIXME -> number of levels nlvdim, and not 1 (done)
! - wind speed and air temp is not available -> introduce (wmeteo)
!
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!
!***********************************************
!---------------------------------------------------------------------
        module toxical
!---------------------------------------------------------------------
        contains
!---------------------------------------------------------------------

	subroutine atoxi3d(it,dt)

! toxi module ARPAV

	use diffusion
	use levels
	use basin
        use shympi
        use elems_dealing
        use ts
        use check
        use conz_util
        use para
        use utility
        use bnd_admin
        use hydro_print
        use initialize
        use bnd_scalar
        use concentration
        use time_util
        use time_admin

	implicit none

	include 'param.h'
	include 'mkonst.h'

	integer it	!time in seconds
	double precision dt		!time step in seconds

	integer nstate
	parameter( nstate = 1 )

	double precision e(nlvdim,nkndim,nstate)	!state vector
	double precision eb(nlvdim,nkndim,nstate)	!boundary vector of state vectors
	save e,eb

        double precision tstot(nstate)              !for mass test

        character*10 what

	integer k,i,l,lmax
	integer itoxi
        integer id
	integer nintp,ivar
	integer nbc
	double precision t,s
	double precision u,v

	integer, save, allocatable :: idtoxi(:)

	double precision eaux(nstate)
	double precision einit(nstate)
	double precision ebound(nstate)
	save einit,ebound

	integer icall,iunit
	integer j
	double precision rlux,rluxaux,itot,fday
	double precision dtt,dttday
	double precision area,vol
	double precision oxysat
	integer ieint,ipint
	integer itanf,nvar
	double precision dtime0,dtime

        integer mode
        double precision ai,lsurf

	logical bcheck
	logical bsedim
	integer ie,ii
	integer kspec
	double precision d
	double precision cbod,nh3,krear,sod
	double precision vel
	double precision windspeed,tempair
	double precision tday,t0,tsec
	double precision stp
        double precision mass
	double precision wsink

	integer iespecial,inspecial
	save iespecial,inspecial
	double precision rkpar,difmol
	save rkpar,difmol
	integer iub,itmcon,idtcon
	save iub,itmcon,idtcon
	integer iubs,itmcons,idtcons
	save iubs,itmcons,idtcons
        save icall

!------------------------------------------------------------------
!	initial and boundary conditions  [mg/l]			??
!
! 	 data einit /0.0, 0., 0.0, 0.0, 0.,   0.,0.,0.0,0.0/
 	 data einit /0.0/
 	 data ebound /0.0/
!
!------------------------------------------------------------------
	data icall /0/
!------------------------------------------------------------------

        what = 'toxi'

!-------------------------------------------------------------------
! initialization
!-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  itoxi = iround(getpar('itoxi'))
	  if( itoxi .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

!         --------------------------------------------------
!	  initialize state variables with einit
!         --------------------------------------------------

	  do k=1,nkn		!loop on nodes
            lmax = ilhkv(k)
            do l=1,lmax
	      do i=1,nstate
	        e(l,k,i) = einit(i)
              end do
	    end do
          end do

!         --------------------------------------------------
!	  initialize state variables from external file
!         --------------------------------------------------

          call inicfil('toxi',e,nstate)

!         --------------------------------------------------
!	  set boundary conditions for all state variables
!         --------------------------------------------------

          nbc = nbnds()
          allocate(idtoxi(nbc))
          idtoxi = 0

	  call get_first_time(itanf)
	  dtime0 = itanf
	  nintp = 2
	  nvar = nstate
          if(bmpi) then
            call bnds_init_mpi(what,dtime0,nintp,nvar,nkn,nlvmax,ebound,idtoxi)
          else
            call bnds_init_new(what,dtime0,nintp,nvar,nkn,nlv,ebound,idtoxi)
          end if

	  !call bnds_init(what,tox3dn,nintp,nstate,nb3dim,toxiarr,ebound)
	  !call bnds_set_def(what,nb3dim,toxiarr) !nvar != 1
	  !call bnds_print('init of '//what,nb3dim,toxiarr)

!         --------------------------------------------------
!	  initialize eco model
!         --------------------------------------------------

!	  call atoxi_ini

!         --------------------------------------------------
!	  parameters for transport/diffusion resolution
!         --------------------------------------------------

          rkpar=getpar('chpar')
          difmol=getpar('difmol')

!         --------------------------------------------------
!	  initialize output 
!         --------------------------------------------------

	  iub = 55
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

          call confop(iub,itmcon,idtcon,nlv,nstate,'tox')

	  write(6,*) 'toxi model initialized...'

	  iubs = 56
          itmcons = iround(getpar('itmcon'))
          idtcons = iround(getpar('idtcon'))

	end if

!-------------------------------------------------------------------
! normal call
!-------------------------------------------------------------------

	kspec = -100
	!kspec = 930
	bcheck = .true.
	bcheck = .false.
	wsink = 0.

!	-------------------------------------------------------------------
!	time management
!	-------------------------------------------------------------------

	t0 = 0.
	call get_timestep(dt)
	tsec = it

!	-------------------------------------------------------------------
!	loop on nodes for biological reactor
!	-------------------------------------------------------------------

        mode = +1               !new time level for volume and depth

	do k=1,nkn		!loop on nodes

          lmax = ilhkv(k)
          !call getmeteo(k,tempair,windspeed)    !meteo FIXME
          !call wmeteo(tempair,windspeed)      !meteo FIXME
          rlux = 1.

          do l=1,lmax
            call dvanode(l,k,mode,d,vol,area)   !gets depth, volume and area
            call getts(l,k,t,s)                 !gets temp and salt
            call getuv(l,k,u,v)                 !gets velocities u/v
            vel = sqrt(u*u+v*v)

            id = 1000*k+l

	    do i=1,nstate
	      eaux(i) = e(l,k,i)
	    end do

	    !call atoxi(id,tsec,dt,d,t,eaux)

	    do i=1,nstate
	      e(l,k,i) = eaux(i)
	    end do
          end do

	end do

!	-------------------------------------------------------------------
!	advection and diffusion
!	-------------------------------------------------------------------

	if( bcheck ) call check_toxi('before advection',e)

	dtime = it
	call bnds_read_new(what,idtoxi,dtime)

!$OMP PARALLEL PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)

	do i=1,nstate

          call scal_adv(what,i,e(1,1,i),idtoxi,rkpar,wsink,difhv,difv,difmol)

          call tsmass (e(1,1,i),1,nlvdim,tstot(i)) !mass control

	end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

	if( bcheck ) call check_toxi('after advection',e)
	write(86,*) it,tstot(1)

!	-------------------------------------------------------------------
!	write of results (file BIO)
!	-------------------------------------------------------------------

	do i=1,nstate
          call confil(iub,itmcon,idtcon,120+i,nlvdim,e(1,1,i))
	end do

!	call toxi_av_shell(e)		!aver/min/max of state vars

	if( bcheck ) call check_toxi('at end',e)

!	-------------------------------------------------------------------
!	debug output
!	-------------------------------------------------------------------

!        k = 100
!        l = 1
!        call getts(l,k,t,s)
!        call writeet(95,it,k,l,e,t,s,nlvdim,nkndim,nstate)

!	-------------------------------------------------------------------
!	end of routine
!	-------------------------------------------------------------------

	end

!*************************************************************

	subroutine writeet(iunit,it,k,l,e,t,s,nlvdim,nknddi,nstate)

! formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nknddi,nstate
	double precision e(nlvdim,nknddi,nstate)
	double precision t
	double precision s

	integer i

	write(iunit,'(i10,11f12.4)') it,(e(l,k,i),i=1,nstate),t,s

	end

!*************************************************************

	subroutine toxi_av_shell(e)

! computes and writes average/min/max of bio variables
!
! id = 260
!
! e(1) average	== 261
! e(1) min	== 262
! e(1) max	== 263
! e(2) average	== 264
! ...
        use para
        use residual

	implicit none

! parameter

	include 'param.h'

	integer nstate
	parameter( nstate = 1 )

	double precision e(nlvdim,nkndim,nstate)	!state vector

! local
	integer idtc,itmc,itsmed
	integer id,nvar
! save
	double precision bioacu(nlvdim,nkndim,nstate)
	double precision biomin(nlvdim,nkndim,nstate)
	double precision biomax(nlvdim,nkndim,nstate)

	integer ivect(8)

	save bioacu,biomin,biomax
	save ivect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nstate

	  id = 260
	  call cmed_init('bav',id,nvar,nlvdim,idtc,itmc,bioacu,biomin,biomax,ivect)

	  icall = 1
	end if

	call cmed_accum(nlvdim,e,bioacu,biomin,biomax,ivect)

	end

!*************************************************************

	subroutine check_toxi(title,e)

! checks bio vars

	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw
        use chk_NaN

	implicit none

	include 'param.h'

	integer nstate
	parameter( nstate = 1 )

	character*(*) title
	double precision e(nlvdim,nkndim,nstate)	!state vector


        character*20 text
	integer i

	write(6,*) 'check_toxi: ',title

        text = '*** bio check e     '
	do i=1,nstate
          write(text(18:19),'(i2)') i
          call check2Dr(nlvdim,nlv,nkn,e(1,1,i),0.d0,1.d+20,text,title)
	end do

	end

!*************************************************************

!---------------------------------------------------------------------
        end module toxical
!---------------------------------------------------------------------
