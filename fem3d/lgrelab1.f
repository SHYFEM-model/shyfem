
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2020  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c revision log :
c
c 20.07.2018	ccf	from scratch
c 25.10.2018	ggu	changed VERS_7_5_51
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 29.01.2020	ggu	old call to nc_output_record_reg() substituted
c 26.06.2021	ggu	better output to terminal
c
c**************************************************************

	subroutine lgrelab

	use clo
	use elabutil
	use elabtime
        use shyfile
        use shyutil
        use shyelab_out
	use shyfem_strings

        use basin
        use mod_depth
        use evgeom
        use levels
        use shympi

	implicit none

	!active particles
        integer, allocatable            :: ida(:)
        integer, allocatable            :: tya(:)
        double precision, allocatable   :: tta(:)
        real, allocatable               :: sa(:)
        integer, allocatable            :: iea(:)
        real, allocatable               :: xa(:),ya(:),za(:)
        integer, allocatable            :: lba(:)
        real, allocatable               :: hla(:)
        real, allocatable               :: ca(:,:)
        real, allocatable               :: agea(:)

	!new inserted particles
        integer, allocatable            :: idn(:)
        integer, allocatable            :: tyn(:)
        double precision, allocatable   :: ttn(:)
        real, allocatable               :: sn(:)
        integer, allocatable            :: ien(:)
        real, allocatable               :: xn(:),yn(:),zn(:)
        integer, allocatable            :: lbn(:)
        real, allocatable               :: hln(:)
        real, allocatable               :: cn(:,:)

	!exited particles
        integer, allocatable            :: ide(:)
        integer, allocatable            :: tye(:)
        double precision, allocatable   :: tte(:)
        real, allocatable               :: se(:)
        integer, allocatable            :: iee(:)
        real, allocatable               :: xe(:),ye(:),ze(:)
        integer, allocatable            :: lbe(:)
        real, allocatable               :: hle(:)
        real, allocatable               :: ce(:,:)

	!mean particles
        integer, allocatable            :: idm(:)
        integer, allocatable            :: tym(:)
        double precision, allocatable   :: ttm(:)
        real, allocatable               :: sm(:)
        integer, allocatable            :: iem(:)
        real, allocatable               :: xm(:),ym(:),zm(:)
        integer, allocatable            :: lbm(:)
        real, allocatable               :: hlm(:)
        real, allocatable               :: cm(:,:)
        real, allocatable               :: agem(:)

        integer, allocatable            :: iaux(:)
        integer, save, allocatable      :: idstore(:)
        real, allocatable               :: paux(:)
        real, save, allocatable         :: ttstore(:)
        integer, allocatable		:: ivars(:)

	integer				:: ncust
	integer				:: date,time
        integer                         :: iwhat
        integer                         :: ierr
	integer				:: n_typ
	integer				:: iu
        integer 			:: i,nvers,lmax,l,idx
        integer 			:: id,idout
        logical                         :: bsphe
        integer                         :: isphe
        integer ifile,ftype,npr,nvar

        INTERFACE
        subroutine lgr_alloc(nn,nc
     +          ,id,ty,tt,s,ie,x,y,z,lb,hl,c)
        integer nn,nc
        integer, allocatable            :: id(:)
        integer, allocatable            :: ty(:)
        double precision, allocatable   :: tt(:)
        real, allocatable               :: s(:)
        integer, allocatable            :: ie(:)
        real, allocatable               :: x(:),y(:),z(:)
        integer, allocatable            :: lb(:)
        real, allocatable               :: hl(:)
        real, allocatable               :: c(:,:)
        end subroutine
        END INTERFACE

        double precision atfirst,atlast,dtime,atime0
        double precision atime,atstart,atnew,atold,atwrite
        character*20 aline
        character*80 header

        integer n_act,n_new,n_ext,n_bea,n_init,n_old
        integer nread,nelab,nrec,nout,nin,nn
	integer iapini,ifileo
	integer ifem_open_file

	header = ' nrec         date_and_time  n_active'
     +                 // '     n_new    n_exit   n_beach    n_type'

c--------------------------------------------------------------

        nread  = 0
        nelab  = 0
	nrec   = 0
        n_init = 0
        ifile  = 0
        id     = 0

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('LGR','lgrelab')

        if( blg2d .and. .not. blgdens ) then
          write(6,*) '-lg2d option only in combination with -lgdens'
          stop 'error stop lgrelab: -lg2d output'
        end if

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------
        call open_new_file(ifile,id,atstart)    !atstart=-1 if no new file
        if( atstart /= -1 ) then
          call dts_format_abs_time(atstart,aline)
          if( .not. bsilent ) then
            write(6,*) 'initial date for next file: ',aline
          end if
        end if

        !--------------------------------------------------------------
        ! set up params and modules
        !--------------------------------------------------------------
        call shy_get_params(id,nkn,nel,npr,nlv,nvar)

        if( .not. bquiet ) call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)
        call mod_depth_init(nkn,nel)
        call basin_set_read_basin(.true.)
        call shy_copy_basin_from_shy(id)
        call shy_copy_levels_from_shy(id)

        call shympi_init(.false.)               !call after basin has been read
        call shympi_set_hlv(nlv,hlv)

        call ev_set_verbose(.not.bquiet)
        call ev_init(nel)
        call set_ev
        call get_coords_ev(isphe)
        bsphe = isphe .eq. 1

        !--------------------------------------------------------------
        ! read number of custom properties
        !--------------------------------------------------------------
        call shy_get_iunit(id,iu)
        read(iu) ncust
	write(*,*)'ncust: ', ncust

        !--------------------------------------------------------------
        ! set up aux arrays, sigma info and depth values
        !--------------------------------------------------------------
        call shyutil_init(nkn,nel,nlv)
        call init_sigma_info(nlv,hlv)
        call shy_make_area
        call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------
        call shy_get_date(id,date,time)
        call dts_to_abs_time(date,time,atime0)
	call elabtime_date_and_time(date,time)	!we work with absolute time
	call elabtime_set_minmax(stmin,stmax)

	call lgr_peek_block_header(iu,dtime,nn,iwhat,ierr)
        call dts_convert_to_atime(datetime_elab,dtime,atime)
	if( ierr /= 0 ) goto 91
        atfirst = atime
        atold = atime

	!--------------------------------------------------------------
	! initilize output for concentration and age
	!--------------------------------------------------------------
        if ( blgdens ) then
          boutput = blgdens
          b2d = blg2d
          nvar = 2
          allocate(ivars(nvar))
          ivars(1) = 84		!custom, density in this case
          ivars(2) = 81		!age
          ftype = 2
          call shyelab_init_output(id,idout,ftype,nvar,ivars)
        end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

        atime = atold
        atwrite = -1

	if( .not. bquiet ) then
          write(6,*)
	  write(6,*)' Particle distrubution over time'
	  write(6,*)
	  write(6,'(a)') header
        end if

	do

          atold = atime

          !----------------------------------------------------------------
          ! read lgr data block --> active particles (plot only these particles)
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,dtime,n_act,iwhat,ierr)
          if( ierr > 0 ) write(6,*) 'error in reading file : ',ierr
          if( ierr /= 0 ) exit
          call lgr_alloc(n_act,ncust
     +          ,ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)
          call lgr_get_block(iu,n_act,ncust,
     +                  ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)


	  call count_lbeach(n_act,tya,n_bea)

          !----------------------------------------------------------------
          ! skip lgr data block --> inserted particles
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,dtime,n_new,iwhat,ierr)
          if( ierr > 0 ) write(6,*) 'error in reading file : ',ierr
          if( ierr /= 0 ) exit
          call lgr_alloc(n_new,ncust
     +          ,idn,tyn,ttn,sn,ien,xn,yn,zn,lbn,hln,cn)
          call lgr_get_block(iu,n_new,ncust,
     +                  idn,tyn,ttn,sn,ien,xn,yn,zn,lbn,hln,cn)

          !----------------------------------------------------------------
          ! skip lgr data block --> exited particles
          !----------------------------------------------------------------
          call lgr_peek_block_header(iu,dtime,n_ext,iwhat,ierr)
          if( ierr > 0 ) write(6,*) 'error in reading file : ',ierr
          if( ierr /= 0 ) exit
          call lgr_alloc(n_ext,ncust
     +          ,ide,tye,tte,se,iee,xe,ye,ze,lbe,hle,ce)
          call lgr_get_block(iu,n_ext,ncust,
     +                  ide,tye,tte,se,iee,xe,ye,ze,lbe,hle,ce)

          !----------------------------------------------------------------
          ! time management
          !----------------------------------------------------------------
	  nrec = nrec + 1
          nread = nread + 1
          call dts_convert_to_atime(datetime_elab,dtime,atime)
          atlast = atime
          if( ierr .ne. 0 ) atnew = atime
          if( elabtime_over_time(atime,atnew,atold) ) exit
          if( .not. elabtime_in_time(atime,atnew,atold) ) cycle
	  !if( .not. elabtime_check_time(atime,atnew,atold) ) cycle

	  call dts_format_abs_time(atime,aline)
	  if( bverb .and. atwrite /= atime ) then
	    write(6,*) 'time : ',atime,'  ',aline
	    atwrite = atime
	  end if

          !----------------------------------------------------------------
          ! Store information and compute particle age
          !----------------------------------------------------------------
          allocate(agea(n_act))
          n_old = n_init
          n_init = n_init + n_new
	  if ( n_new > 0 ) then
            if( n_init == n_new ) then
               allocate(ttstore(n_init))
               allocate(idstore(n_init))
               ttstore = ttn
               idstore = idn
            else
               allocate(paux(n_init))
               allocate(iaux(n_init))
               paux(1:n_old) = ttstore(1:n_old)
               call move_alloc(paux,ttstore)
               ttstore(n_old+1:n_init) = ttn(1:n_new)
               iaux(1:n_old) = idstore(1:n_old)
               call move_alloc(iaux,idstore)
               idstore(n_old+1:n_init) = idn(1:n_new)
            end if
          end if

          do i = 1,n_act
            idx = minloc(abs(idstore - ida(i)), 1)
            agea(i) = tta(i) - ttstore(idx)
          end do
          agea = max(agea,0.)
          agea = agea / 86400.		!from sec to day

          !----------------------------------------------------------------
          ! print lgr informations
          !----------------------------------------------------------------

          n_typ = 0
          if ( n_act > 0 ) n_typ = maxval(tya)

	  if( .not. bquiet ) 
     +         write(6,100) nrec,aline,n_act,n_new,n_ext,n_bea,n_typ

          !----------------------------------------------------------------
          ! Compute mean position for active particle based on type
          !----------------------------------------------------------------
          if ( blgmean ) then
            if (n_typ > 0 ) then
              allocate(agem(n_typ))
              call lgr_alloc(n_typ,ncust
     +          ,idm,tym,ttm,sm,iem,xm,ym,zm,lbm,hlm,cm)
              call lgr_mean_posit(n_act,ncust,tya,agea,sa,xa,ya,hla,ca,
     +                       n_typ,tym,agem,sm,xm,ym,hlm,cm)
              ttm = atime 

              call lgr_write_mean_posit(n_typ,ncust,aline,
     +                       agem,sm,xm,ym,hlm,cm)

              deallocate(agem)
            end if
          end if

          !----------------------------------------------------------------
          ! Compute active particle density on nodes and write .lgc.shy
          !----------------------------------------------------------------
          if ( blgdens ) then
            dtime = atime - atime0
            call lgr_output_conc(idout,ftype,nvar,ivars,dtime,
     +                  bsphe,n_act,iea,xa,ya,hla,lba,tya,agea)
          end if

          deallocate(agea)
	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

	  write(6,'(a)') header

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*)
	  call dts_format_abs_time(atfirst,aline)
          write(6,*) 'first time record: ',aline
	  call dts_format_abs_time(atlast,aline)
          write(6,*) 'last time record:  ',aline

	  write(6,*)
	  write(6,*) nread,' data records read'
	  write(6,*) nrec ,' time records read'
	  write(6,*) nelab,' time records elaborated'
	  write(6,*)
	end if

        call shyelab_final_output(id,idout,nvar,ivars)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	return
   85	continue
	write(6,*) 'it,itvar,i '
	stop 'error stop extelab: time mismatch'
   91	continue
	write(6,*) 'error reading first data record'
	write(6,*) 'maybe the file is empty'
	stop 'error stop extelab: empty record'
   93	continue
	write(6,*) 'error reading header of file'
	write(6,*) 'maybe the file is empty'
	stop 'error stop extelab: empty header'
   99	continue
	write(6,*) 'error writing to file unit: '
	stop 'error stop extelab: write error'
  100   format(i5,a22,5i10)
	end

c***************************************************************

        subroutine count_lbeach(n,ty,nb)

	implicit none

	integer, intent(in)	:: n
	integer, intent(in)	:: ty(n)
	integer, intent(out)	:: nb

	integer			:: i

	nb = 0
	do i = 1,n
	  if (ty(i) < 0 ) nb = nb + 1	
	end do

	end subroutine count_lbeach
	
!***************************************************************
! Write mean particles position and characteristics

        subroutine lgr_write_mean_posit(nm,nc,aline,
     +                       agem,sm,xm,ym,hlm,cm)

        implicit none

        integer, intent(in)           :: nm             !number of types
        integer, intent(in)           :: nc 
        character*20 		      :: aline
        real, intent(in)              :: agem(nm)
        real, intent(in)              :: sm(nm)
        real, intent(in)              :: xm(nm),ym(nm)
        real, intent(in)              :: hlm(nm)
        real, intent(in)              :: cm(nm,nc)

	integer			      :: i,j,iu
        integer, allocatable,save     :: ium(:)
        integer, allocatable          :: iaux(:)
        character*80		      :: name,format
        character*5		      :: numb
        character*18		      :: basen='lagrange_mean_traj'

        integer, save     	      :: icall = 0
        integer, save                 :: nty = 0

!       ------------------------------------------------------
!       initialization of output files
!       ------------------------------------------------------
        if( icall .eq. 0 ) then
          allocate(ium(nm))
          do i = 1,nm
            write(numb,'(i5)') i
            numb = adjustl(numb)
            name = trim(basen)// '.' //numb
            call get_new_unit(iu)
            open(iu,file=name,form='formatted',status='unknown')
            ium(i) = iu
          end do
          nty = nm
          icall = 1
        end if

!       ------------------------------------------------------
!       normal call
!       ------------------------------------------------------
        if ( nm > nty ) then
          allocate(iaux(nm))
          iaux(1:nty) = ium(1:nty)
          call move_alloc(iaux,ium)
          do i = nty+1,nm
            write(numb,'(i5)') i
            numb = adjustl(numb)
            name = trim(basen)// '.' //numb
            call get_new_unit(iu)
            open(iu,file=name,form='formatted',status='unknown')
            ium(i) = iu
          end do
          nty = nm
        end if

!       --------------------------------------------------
!       Write to output files
!       --------------------------------------------------
        write(format,'(a,i3,a)') '(a20,2f14.5,3f10.3,',nc,'f10.3)'

	do i = 1,nm
           iu = ium(i)
           write(iu,format) aline,xm(i),ym(i),hlm(i),agem(i),sm(i),
     +                   (cm(i,j),j=1,nc)
        end do

        end subroutine lgr_write_mean_posit

!***************************************************************
! writes a .lgc.shy file containing:
!   - particles as density (concentration) on nodes
!   - mean age of particles on nodes
! if the reg option is used, then the density and age calculations
! are performed on regular grid

        subroutine lgr_output_conc(idout,ftype,nvar,ivars,dtime,
     +			bsphe,na,iea,x,y,hl,lb,ty,age)

        use evgeom
        use basin
        use levels
        use elabutil
        use shyelab_out

        implicit none

        integer, intent(in)		:: idout
        integer, intent(in)		:: ftype
        integer, intent(in)		:: nvar
        integer, intent(in)		:: ivars(nvar)
        double precision, intent(in)	:: dtime
        logical, intent(in)		:: bsphe
        integer, intent(in)             :: na             !number of particles
        integer, intent(in)             :: iea(na)
        real, intent(in)                :: x(na)
        real, intent(in)                :: y(na)
        real, intent(in)                :: hl(na)
        integer, intent(in)             :: lb(na)
        integer, intent(in)             :: ty(na)
        real, intent(in)                :: age(na)

        integer, allocatable		:: ecount(:,:)
        integer, allocatable		:: kcount(:,:)
        real, allocatable		:: area(:,:)
        real, allocatable		:: density(:,:)
        real, allocatable		:: eage(:,:)
        real, allocatable		:: kage(:,:)

        real, parameter :: pi=3.14159265358979323846, rad = pi/180.
	real				:: xp,x1,x2,x3,x4
	real				:: yp,y1,y2,y3,y4
	real				:: dx,dy 
	integer				:: n,m,ix,iy
        character*60 string
	integer iunit,ncid,nvers,var_id
        integer ie,ii,k,l,lmax,np,ll
        integer ic,i
        real area_node,aage
        real rlat,flon,flat

!---------------------------------------------------------
! initialize
!---------------------------------------------------------

	if ( blg2d ) then
          lmax = 1
	  nlvdi = 1
	else 
          lmax = nlv
	end if

        if ( breg ) then
          np = nxreg * nyreg
        else
          np = nkn
        end if

        allocate(ecount(nlvdi,nel))
        allocate(eage(nlvdi,nel))
        allocate(kcount(nlvdi,np))
        allocate(area(nlvdi,np))
        allocate(density(nlvdi,np))
        allocate(kage(nlvdi,np))

        ecount = 0
        area   = 0.
        density = 0.
        eage   = 0.
        kcount = 0
        kage   = 0.
        flon   = 0.
        flat   = 0.

!---------------------------------------------------------
! regular --> count particles on regular grid
!---------------------------------------------------------

	if ( breg ) then

          n = nkn          
          m = nyreg          
          do i=1,na
            ie = iea(i)
            if ( ie < 1 ) cycle
            l  = lb(i)
            if ( l < 1 ) l = ilhv(ie)	!count particle on bottom
            if ( blg2d ) l = 1
            xp = x(i)
            yp = y(i)
            ix = minloc(abs(xlon-xp),1)
            iy = minloc(abs(ylat-yp),1)
            ii = ix + (iy-1)*nxreg
            kcount(l,ii) = kcount(l,ii) + 1
            kage(l,ii) = kage(l,ii) + age(i)
          end do

          ii = 0
          do iy=1,nyreg
            if ( bsphe ) then
              rlat = ylat(iy)*rad
              flon=111412.84*cos(rlat) - 93.5*cos(3.*rlat)
              flat=111132.92 - 559.82*cos(2.*rlat) +1.175*cos(4.*rlat)
            endif
            dx = regpar(5)*flon
            dy = regpar(6)*flat
            area_node = dx*dy
            do ix=1,nxreg
              ii = ii + 1
              if ( hcoord(ii) == flag ) then
                density(:,ii) = flag
                kage(:,ii) = flag
              else
                do l = 1,nlvdi
                  if( kcount(l,ii) > 0 ) then
                    density(l,ii) =  kcount(l,ii) / area_node
                    kage(l,ii) = kage(l,ii) / kcount(l,ii)
                  end if
                end do
              end if
            end do
          end do

	else

!---------------------------------------------------------
! irregular --> count particles in elements and convert to nodes
!---------------------------------------------------------
          n = nkn          
          m = 1          

          do i=1,na
            ie = iea(i)
            if ( ie < 1 ) cycle
            l = lb(i)
	    if ( l < 1 ) l = ilhv(ie)
            if ( blg2d ) l = 1
            ecount(l,ie) = ecount(l,ie) + 1
            eage(l,ie) = eage(l,ie) + age(i)
          end do

          do ie=1,nel
            ll = min(lmax,ilhv(ie))
  	    do l = 1,ll
              ic = ecount(l,ie)
              if ( ic > 0 ) then
                area_node = 4*ev(10,ie)
                aage = eage(l,ie)/ic
                do ii=1,3
                  k = nen3v(ii,ie)
                  kcount(l,k) = kcount(l,k) + ic
                  kage(l,k) = kage(l,k) + aage*area_node
                  area(l,k) = area(l,k) + area_node
                end do
              end if
            end do
          end do

          do k=1,nkn
            ll = min(lmax,ilhkv(k))
            do l = 1,ll
              area_node = area(l,k)
              if ( kcount(l,k) > 0 ) then
                density(l,k) = kcount(l,k) / area_node
                kage(l,k) = kage(l,k) / area_node
              end if
            end do
          end do
 
        end if
	 
c---------------------------------------------------------
c write to file
c---------------------------------------------------------
        call shyelab_header_output(idout,ftype,dtime,nvar)

        if( outformat == 'shy' .or. outformat == 'native' ) then
          call shy_write_output_record(idout,dtime,ivars(1),n,m
     +                                        ,lmax,nlvdi,density)
          call shy_write_output_record(idout,dtime,ivars(2),n,m
     +                                        ,lmax,nlvdi,kage)
        else if( outformat == 'gis' ) then
          call gis_write_record(dtime,ivars(1),np,lmax,ilcoord
     +                                ,density,xcoord,ycoord)
          call gis_write_record(dtime,ivars(2),np,lmax,ilcoord
     +                                ,kage,xcoord,ycoord)
        else if( outformat == 'fem' ) then
          iunit = idout
          nvers = 0
          call get_string_description(ivars(1),string)
          call fem_file_write_data(iformat,iunit
     +                        ,nvers,np,lmax
     +                        ,string
     +                        ,ilcoord,hcoord
     +                        ,nlvdi,density)
          call get_string_description(ivars(2),string)
          call fem_file_write_data(iformat,iunit
     +                        ,nvers,np,lmax
     +                        ,string
     +                        ,ilcoord,hcoord
     +                        ,nlvdi,kage)
        else if( outformat == 'nc' ) then
          ncid = idout
          var_id = var_ids(1)
          call nc_output_record(ncid,var_id,np,density)
          var_id = var_ids(2)
          call nc_output_record(ncid,var_id,np,kage)
        else if( outformat == 'off' ) then
          ! nothing to be done
        else
          write(6,*) 'outformat = ',trim(outformat)
          stop 'error stop: outformat not recognized'
        end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end subroutine lgr_output_conc

!***************************************************************

