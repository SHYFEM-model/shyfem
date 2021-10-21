
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-1999,2001,2003,2007-2008,2010  Georg Umgiesser
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2018  Petras Zemlys
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

! revision log :
!
! 18.11.1998	ggu	check dimensions with dimnos
! 06.04.1999	ggu	some cosmetic changes
! 03.12.2001	ggu	some extra output -> place of min/max
! 09.12.2003	ggu	check for NaN introduced
! 07.03.2007	ggu	easier call
! 08.11.2008	ggu	do not compute min/max in non-existing layers
! 07.12.2010	ggu	write statistics on depth distribution (depth_stats)
! 06.05.2015	ggu	noselab started
! 05.06.2015	ggu	many more features added
! 10.09.2015	ggu	std and rms for averaging implemented
! 11.09.2015	ggu	write in gis format
! 23.09.2015	ggu	handle more than one file (look for itstart)
! 16.10.2015	ggu	started shyelab
! 25.05.2016	ggu	changed VERS_7_5_10
! 30.05.2016	ggu	changed VERS_7_5_11
! 10.06.2016	ggu	shyplot now plots fem files
! 13.06.2016	ggu	shyplot now plots barotropic vars (layer==0)
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.10.2016	ggu	changed VERS_7_5_19
! 31.10.2016	ggu	shyplot restructured... directional plot still broken
! 12.01.2017	ggu	changed VERS_7_5_21
! 20.01.2017	ggu	changed VERS_7_5_22
! 13.02.2017	ggu	changed VERS_7_5_23
! 14.02.2017	ggu	bug fix in plotting regular fem files - introduced il
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 11.07.2017	ggu	changed VERS_7_5_30
! 09.10.2017	ggu	changed VERS_7_5_33
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	shyplot unified and simplified for output
! 17.11.2017	ggu	changed VERS_7_5_37
! 17.11.2017	ggu	changed VERS_7_5_38
! 24.01.2018	ggu	changed VERS_7_5_41
! 07.06.2018	pzy	new module plot_fonts for font size definition
! 21.06.2018	ccf	shyplot working also for lagrangian particles
! 06.07.2018	ggu	changed VERS_7_5_48
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 21.05.2019	ggu	changed VERS_7_5_62
! 08.06.2021	ggu	forgot to call populate_strings()
! 25.06.2021	ggu	populate_strings() before plotutil_init()
! 25.06.2021	ggu	in plot_lgr_file() call shympi_init() after basin init
! 21.10.2021	ggu	fixed for vertical velocity as overlay
!
! notes :
!
! for customization please see file supcust.f
!
!**************************************************************

	program shyplot

	use plotutil

	implicit none

	call shyfem_copyright('shyplot - plotting SHYFEM files')

	call populate_strings

	call plotutil_init('SHY')
	call classify_files

	if( lgrfilename /= ' ' ) then
	  call plot_lgr_file
	else if( shyfilename /= ' ' ) then
	  call plot_shy_file
	else if(femfilename /= ' ' ) then
	  call plot_fem_file
	else if( basfilename /= ' ' ) then
	  call plot_bas_file
	end if

	end

!**************************************************************

	subroutine plot_bas_file

        use mod_depth
        use mod_hydro_plot
        use mod_geom
        use evgeom
        use basin
        use plotutil
        use shympi

	implicit none

	integer ivar

        call read_command_line_file(basfilename)

	call shympi_init(.false.)

	call bash_verbose(bsdebug)
	call ev_set_verbose(.not.bquiet)
        call ev_init(nel)
        call set_ev

        call mod_geom_init(nkn,nel,ngr)
        call set_geom

        call mod_depth_init(nkn,nel)
        call mod_hydro_plot_init(nkn,nel,1,nel)

        call makehev(hev)
        call makehkv(hkv)
        call allocate_2d_arrays

        call init_plot

	ivar = 5			!bathymetry
	call init_nls_fnm
	call read_str_files(-1)
	call read_str_files(ivar)
        call initialize_color

	if( .not. bquiet ) write(6,*) 'plotting basin...'
        call qopen
	call plobas
        call qclose

	end

!**************************************************************

	subroutine plot_lgr_file

        use clo
        !use elabutil
        use plotutil
        use elabtime
        use shyfile
        use shyutil
        use shympi

        use basin
        use levels
        use evgeom
        use mod_depth
        use mod_geom
        use mod_hydro_plot

	implicit none

        integer iunit,i,n,nvers,lmax,l
        integer, allocatable            :: ida(:)
        integer, allocatable            :: tya(:)
        double precision, allocatable   :: tta(:)
        real, allocatable               :: sa(:)
        integer, allocatable            :: iea(:)
        real, allocatable               :: xa(:),ya(:),za(:)
        integer, allocatable            :: lba(:)
        real, allocatable               :: hla(:)
        real, allocatable               :: ca(:,:)
        integer, allocatable            :: aplot(:)
        real, allocatable               :: ra(:)
        real, allocatable               :: agea(:)

        integer, allocatable            :: idn(:)
        integer, allocatable            :: tyn(:)
        double precision, allocatable   :: ttn(:)
        real, allocatable               :: sn(:)
        integer, allocatable            :: ien(:)
        real, allocatable               :: xn(:),yn(:),zn(:)
        integer, allocatable            :: lbn(:)
        real, allocatable               :: hln(:)
        real, allocatable               :: cn(:,:)

        integer, allocatable            :: idm(:)
        integer, allocatable            :: tym(:)
        double precision, allocatable   :: ttm(:)
        real, allocatable               :: sm(:)
        integer, allocatable            :: iem(:)
        real, allocatable               :: xm(:),ym(:),zm(:)
        integer, allocatable            :: lbm(:)
        real, allocatable               :: hlm(:)
        real, allocatable               :: cm(:,:)
        integer, allocatable            :: mplot(:)
        real, allocatable               :: rm(:)
        real, allocatable               :: agem(:)

        real, save,  allocatable        :: xall(:,:),yall(:,:),rall(:,:)
        real, save, allocatable         :: xmll(:,:),ymll(:,:),rmll(:,:)
        real, allocatable               :: aux(:,:)

        integer, allocatable 		:: iaux(:)
        integer, save, allocatable 	:: idstore(:)
        real, allocatable 		:: paux(:)
        real, save, allocatable 	:: ttstore(:)

        integer                         :: iwhat
        integer                         :: ierr
	integer				:: lgmean
        logical				:: btraj = .false.
        logical				:: blgmean = .false.

        integer, save 		        :: nn_old,nt_old,na_old
        integer n_act,n_new,n_ext,n_init,n_typ
        integer ncust
        character*80 name
        logical ptime_ok,ptime_end
        integer irec,nplot,idx
        integer nb,nout,np
        integer ifileo
        integer getlev,getvar
        real rmin,rmax,flag
        integer nt
	character*80 line
        logical bhasbasin
	integer isphe
        real getpar
        double precision dgetpar
	integer id,idold,npr,nvar,ftype

	integer isub
	logical bsect,bskip
        integer date,time
        integer datetime(2)
        double precision atime,atime0,dtime

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

        !--------------------------------------------------------------
        ! set command line parameters
        !--------------------------------------------------------------

        call init_nls_fnm
        call read_str_files(-1)
        call read_str_files(ivar3)

        btraj  = nint(getpar('lgrtrj')) == 1
        lgmean  = nint(getpar('lgmean'))
        blgmean = lgmean > 0

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

        id = 0
        idold = 0
	call open_next_file_by_name(lgrfilename,idold,id)
	if( id == 0 ) stop

	!--------------------------------------------------------------
	! set up params and modules
	!--------------------------------------------------------------

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	if( .not. bquiet ) call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)

	call basin_set_read_basin(.true.)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)
	call bash_verbose(bsdebug)

        call shympi_init(.false.)

        call mod_depth_init(nkn,nel)
	call allocate_2d_arrays

        isphe = nint(getpar('isphe'))
        call set_coords_ev(isphe)
	call ev_set_verbose(.not.bquiet)
        call set_ev
        call set_geom
        call get_coords_ev(isphe)
        call putpar('isphe',float(isphe))

	!--------------------------------------------------------------
	! read ncust
	!--------------------------------------------------------------

        call shy_get_iunit(id,iunit)
        read(iunit) ncust

	!--------------------------------------------------------------
	! set time
	!--------------------------------------------------------------

        call ptime_init
	call shy_get_date(id,date,time)
        call dts_to_abs_time(date,time,atime0)
        call ptime_set_date_time(date,time)
        call elabtime_date_and_time(date,time)
        call elabtime_set_minmax(stmin,stmax)
	call elabtime_set_inclusive(.false.)

        irec = 0
	nplot = 0
	n_init = 0
        flag = dflag

        !--------------------------------------------------------------
        ! initialize plot
        !--------------------------------------------------------------

        call initialize_color
        call qopen

        if (btraj) then
	  write(6,*) 'Plotting lagrangian particles trajectories'
        else
	  write(6,*) 'Plotting lagrangian particles positions'
        end if

	!----------------------------------------------------------------
	! set what to plot with color with option varnam in plots
	!----------------------------------------------------------------

        call ivar2string(ivar3,name,isub)
	if( ivar3 == 0 ) name = 'lgr'
	!write(6,*) ivar3,' ',trim(name),isub

        if ( name .eq. 'lgr' ) then           !nothing
           write(6,*)'No color used...'
        else if( name .eq. 'lagtyp' ) then  	!type of particle
           write(6,*)'Variable to be plotted: particle type'
        else if( name .eq. 'lagdep' ) then  	!absolute depth
           write(6,*)'Variable to be plotted: particle depth'
        else if( name .eq. 'lagage' ) then    !age [d]
           write(6,*)'Variable to be plotted: particles age'
        else if( name .eq. 'lagcus' ) then	!custom
           write(6,*)'Variable to be plotted: particle custom prop.'
        else
           goto 99
        end if

	!--------------------------------------------------------------
	! loop on data (time loop)
	!--------------------------------------------------------------

	do

	  !----------------------------------------------------------------
	  ! read lgr data block --> active particles (plot only these particles)
	  !----------------------------------------------------------------

          call lgr_peek_block_header(iunit,dtime,n_act,iwhat,ierr)
          if( ierr /= 0 ) exit
          call lgr_alloc(n_act,ncust
     +          ,ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)
          allocate(aplot(n_act))
          allocate(ra(n_act))
          call lgr_get_block(iunit,n_act,ncust,
     +                  ida,tya,tta,sa,iea,xa,ya,za,lba,hla,ca)
          n = n_act

	  !----------------------------------------------------------------
	  ! skip lgr data block --> inserted particles
	  !----------------------------------------------------------------

          call lgr_peek_block_header(iunit,dtime,n_new,iwhat,ierr)
          if( ierr /= 0 ) exit
	  if ( n_new == 0 ) then
            call lgr_skip_block(iunit,n_new,ncust)
          else

 	    !----------------------------------------------------------------
	    ! read inserted particles and store information for age
	    !----------------------------------------------------------------

            call lgr_alloc(n_new,ncust
     +          ,idn,tyn,ttn,sn,ien,xn,yn,zn,lbn,hln,cn)
            call lgr_get_block(iunit,n_new,ncust,
     +                  idn,tyn,ttn,sn,ien,xn,yn,zn,lbn,hln,cn)
            nn_old = n_init
            n_init = n_init + n_new
            if( n_init == n_new ) then
               allocate(ttstore(n_init))
               allocate(idstore(n_init))
               ttstore = ttn
               idstore = idn
            else
               allocate(paux(n_init))
               allocate(iaux(n_init))
               paux(1:nn_old) = ttstore(1:nn_old)
               call move_alloc(paux,ttstore)
               ttstore(nn_old+1:n_init) = ttn(1:n_new)
               iaux(1:nn_old) = idstore(1:nn_old)
               call move_alloc(iaux,idstore)
               idstore(nn_old+1:n_init) = idn(1:n_new)
            end if
          end if

	  !----------------------------------------------------------------
	  ! skip lgr data block --> exites particles
	  !----------------------------------------------------------------

          call lgr_peek_block_header(iunit,dtime,n_ext,iwhat,ierr)
          if( ierr /= 0 ) exit
          call lgr_skip_block(iunit,n_ext,ncust)

          if (n_act > 0 ) irec = irec + 1

	  !----------------------------------------------------------------
	  ! compute age
	  !----------------------------------------------------------------

          allocate(agea(n_act))
          do i = 1,n_act
            idx = minloc(abs(idstore - ida(i)), 1)
            agea(i) = tta(i) - ttstore(idx)
          end do
          agea = agea / 86400.

          !----------------------------------------------------------------
          ! Compute mean position for active particle based on type
          ! assuming that type starts from 1 
          !----------------------------------------------------------------

          n_typ = 0
          if ( blgmean ) n_typ = maxval(tya)
          allocate(agem(n_typ))
          allocate(mplot(n_typ))
          allocate(rm(n_typ))
          call lgr_alloc(n_typ,ncust
     +        ,idm,tym,ttm,sm,iem,xm,ym,zm,lbm,hlm,cm)

          if ( blgmean ) then
            call lgr_mean_posit(n_act,ncust,tya,agea,sa,xa,ya,hla,ca,
     +                     n_typ,tym,agem,sm,xm,ym,hlm,cm)
            ttm = atime
          end if

	  !----------------------------------------------------------------
	  ! set what to plot with color with option varnam in plots
	  !----------------------------------------------------------------

          if ( name .eq. 'lgr' ) then           !nothing
             ra = 0.
             rm = 0.
          else if( name .eq. 'lagtyp' ) then  	!type of particle
             ra = tya
             rm = tym
          else if( name .eq. 'lagdep' ) then  	!absolute depth [m]
             ra = hla
             rm = hlm
          else if( name .eq. 'lagage' ) then    !age [day]
             ra = agea
             rm = agem
          else if( name .eq. 'lagcus' ) then	!custom
             ra = ca(:,1)
             rm = cm(:,1)
          else
              goto 99
          end if

          rmax = maxval(ra)
          rmin = minval(ra)

          !----------------------------------------------------------------
          ! for plotting the trajectories allocate xall, yall and rall
	  ! CCF ACCOUNT FOR VARIABLE NUMBER OF PARTICLES IN TIME, TO BE DONE
          !----------------------------------------------------------------

          if (btraj) then
            if (irec == 1) then
              nt = 50
              allocate(xall(n,0:nt))
              allocate(yall(n,0:nt))
              allocate(rall(n,0:nt))
              xall(:,0) = xa
              yall(:,0) = ya
              rall(:,0) = ra
              allocate(xmll(n_typ,0:nt))
              allocate(ymll(n_typ,0:nt))
              allocate(rmll(n_typ,0:nt))
              xmll(:,0) = xm
              ymll(:,0) = ym
              rmll(:,0) = rm
              na_old = n
              nt_old = n_typ
            end if

            if ( n > na_old ) then
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:na_old,0:irec-1) = xall(1:na_old,0:irec-1)
              call move_alloc(aux,xall)
              xall(na_old+1:n,0:irec-1) = flag
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:na_old,0:irec-1) = yall(1:na_old,0:irec-1)
              call move_alloc(aux,yall)
              yall(na_old+1:n,0:irec-1) = flag
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:na_old,0:irec-1) = rall(1:na_old,0:irec-1)
              call move_alloc(aux,rall)
              rall(na_old+1:n,0:irec-1) = flag
              na_old = n
            end if
  
            if ( n_typ > nt_old ) then
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:nt_old,0:irec-1) = xmll(1:nt_old,0:irec-1)
              call move_alloc(aux,xmll)
              xmll(nt_old+1:n_typ,0:irec-1) = flag
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:nt_old,0:irec-1) = ymll(1:nt_old,0:irec-1)
              call move_alloc(aux,ymll)
              ymll(nt_old+1:n_typ,0:irec-1) = flag
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:nt_old,0:irec-1) = rmll(1:nt_old,0:irec-1)
              call move_alloc(aux,rmll)
              rmll(nt_old+1:n_typ,0:irec-1) = flag
              nt_old = n_typ
            end if

            if ( irec == nt ) then
              nt = nt + 50
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:n,0:irec-1) = xall(1:n,0:irec-1)
              call move_alloc(aux,xall)
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:n,0:irec-1) = yall(1:n,0:irec-1)
              call move_alloc(aux,yall)
              allocate(aux(n,0:nt))
              aux = 0
              aux(1:n,0:irec-1) = rall(1:n,0:irec-1)
              call move_alloc(aux,rall)
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:n_typ,0:irec-1) = xmll(1:n_typ,0:irec-1)
              call move_alloc(aux,xmll)
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:n_typ,0:irec-1) = ymll(1:n_typ,0:irec-1)
              call move_alloc(aux,ymll)
              allocate(aux(n_typ,0:nt))
              aux = 0
              aux(1:n_typ,0:irec-1) = rmll(1:n_typ,0:irec-1)
              call move_alloc(aux,rmll)
            end if
  
            !----------------------------------------------------------------
            ! store x,y and r
            !----------------------------------------------------------------

            xall(:,irec) = xa
            yall(:,irec) = ya
            rall(:,irec) = ra
            xmll(:,irec) = xm
            ymll(:,irec) = ym
            rmll(:,irec) = rm
          end if

          !----------------------------------------------------------------
          ! set vertical level to plot with option layer in shyplot
          !----------------------------------------------------------------

          mplot = 1
          aplot = 1
          call setlev(layer)
        
          np = n 
          if ( layer /= 0 ) then
            aplot = 0
            np = 0.
            do i = 1,n
              l = lba(i)
              if ( l == layer ) then
                aplot(i) = 1
                np = np + 1
              end if
            end do
          end if

	  !----------------------------------------------------------------
	  ! plot particles if time in timerange
	  !----------------------------------------------------------------

	  atime = dtime + atime0
          call dts_format_abs_time(atime,line)
          call ptime_set_atime(atime)

          bskip = .false.
          if( elabtime_over_time(atime) ) exit
          if( .not. elabtime_in_time(atime) ) bskip = .true.
          if( ifreq > 0 .and. mod(irec,ifreq) /= 0 ) bskip = .true.
          if( n_act == 0 ) bskip = .true.

          if( bskip ) then
            if( bverb ) then
              write(6,*) irec,trim(line),'   ...skipping'
            end if
            deallocate(aplot)
            deallocate(ra)
            deallocate(agea)
            deallocate(agem)
            deallocate(rm)
            deallocate(mplot)
	    cycle
          else
            if( .not. bsilent ) then
              write(6,*) irec,trim(line),'   ...plotting'
            end if
          end if
 
          !call reset_dry_mask

          if( bverb ) write(6,*) 'plotting particles ',atime,np
          if( bverb ) write(6,*) 'plotting rlag: ',trim(name),rmin,rmax

          if (btraj) then
            call plo_traj(n,nt,irec,lgmean,xall,yall,rall,aplot,
     +              n_typ,xmll,ymll,rmll,mplot,'trajectories')
          else
            call plo_part(n,xa,ya,ra,aplot,'particles')
          end if


          nplot = nplot + 1

          deallocate(aplot)
          deallocate(ra)
          deallocate(agea)
          deallocate(agem)
          deallocate(rm)
          deallocate(mplot)

        enddo

!--------------------------------------------------------------
! end of loop on data (time loop)
!--------------------------------------------------------------

        call qclose

        if( .not. bsilent ) then
          write(6,*) 'total number of plots: ',nplot
        end if

        return
   99   continue
        write(6,*) 'Unknown varnam: ',name
        stop 'error stop plot_lgr_file: varnam error'

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!**************************************************************

	subroutine plot_shy_file

	use clo
	!use elabutil
	use plotutil
	use elabtime
	use shyfile
	use shyutil

        use basin
        use levels
        use evgeom
        use mod_depth
	use mod_hydro_plot
        use mod_hydro
        use mod_hydro_vel
        use mod_hydro_print

	implicit none

	real, allocatable :: cv2(:)
	real, allocatable :: cv3(:,:)
	real, allocatable :: cv3all(:,:,:)

	integer, allocatable :: idims(:,:)
	integer, allocatable :: il(:)
	integer, allocatable :: ivars(:)
	character*80, allocatable :: strings(:)

	logical bhydro,bscalar,bsect,bvect,bvel,bcycle
	logical bregplot,bregdata
	logical btime
	integer nx,ny
	integer irec,nplot,nread,nin,nold
	integer nvers
	integer nvar,npr
	integer ierr,ivel,iarrow
	integer it
	integer ivar,iaux,nv
	integer ivarplot(2),ivs(2)
	integer iv,i,j,l,k,lmax,node
	integer ip,np
	integer ifile,ftype
	integer id,idout,idold
	integer n,m,nndim,nn
	integer naccum
	integer isphe
	integer date,time
	character*80 title,name,file
	character*80 basnam,simnam,varline
	real rnull
	real cmin,cmax,cmed,vtot
	real dx,dy
	double precision dtime
	double precision atime

	integer iapini
	integer ifem_open_file
	integer getisec
	real getpar

!--------------------------------------------------------------
! initialize everything
!--------------------------------------------------------------

	irec=0
	nplot=0
	nread=0
	rnull=0.
	rnull=-1.
	bzeta = .false.		!file has zeta information
	ifile = 0
	id = 0
	idold = 0
	bregdata = .false.	!shy file data is not regular

	!--------------------------------------------------------------
	! set command line parameters
	!--------------------------------------------------------------

	call init_nls_fnm
	call read_str_files(-1)
	call read_str_files(ivar3)

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call open_next_file_by_name(shyfilename,idold,id)
	if( id == 0 ) stop

	!--------------------------------------------------------------
	! set up params and modules
	!--------------------------------------------------------------

	call populate_strings

	call shy_get_params(id,nkn,nel,npr,nlv,nvar)
	call shy_get_ftype(id,ftype)

	if( .not. bquiet ) call shy_info(id)

        call basin_init(nkn,nel)
        call levels_init(nkn,nel,nlv)

	call basin_set_read_basin(.true.)
	call shy_copy_basin_from_shy(id)
	call shy_copy_levels_from_shy(id)
	call bash_verbose(bsdebug)

        call mod_depth_init(nkn,nel)
	call allocate_2d_arrays
	np = nel
	call mod_hydro_plot_init(nkn,nel,nlv,np)

	call mod_hydro_init(nkn,nel,nlv)
	call mod_hydro_vel_init(nkn,nel,nlv)
	call mod_hydro_print_init(nkn,nlv)
	znv = 0.
	zenv = 0.

        isphe = nint(getpar('isphe'))
        call set_coords_ev(isphe)
	call ev_set_verbose(.not.bquiet)
        call set_ev
        call set_geom
        call get_coords_ev(isphe)
        call putpar('isphe',float(isphe))

	!--------------------------------------------------------------
	! set time
	!--------------------------------------------------------------

        call ptime_init
	call shy_get_date(id,date,time)
        call ptime_set_date_time(date,time)
        call elabtime_date_and_time(date,time)
        call elabtime_set_minmax(stmin,stmax)
	call elabtime_set_inclusive(.false.)

	!--------------------------------------------------------------
	! set dimensions and allocate arrays
	!--------------------------------------------------------------

	bhydro = ftype == 1
	bscalar = ftype == 2

	if( bhydro ) then		!OUS
	  if( nvar /= 4 ) goto 71
	  nndim = 3*nel
	  allocate(il(nel))
	  il = ilhv
	else if( bscalar ) then		!NOS
	  nndim = nkn
	  allocate(il(nkn))
	  il = ilhkv
	else
	  goto 76	!relax later
	end if

	allocate(cv2(nndim))
	allocate(cv3(nlv,nndim))
	allocate(cv3all(nlv,nndim,0:nvar))
	allocate(idims(4,nvar))
	allocate(ivars(nvar),strings(nvar))

	!--------------------------------------------------------------
	! set up aux arrays, sigma info and depth values
	!--------------------------------------------------------------

	call shyutil_init(nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	call shy_make_area
	call outfile_make_depth(nkn,nel,nen3v,hm3v,hev,hkv)

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	bminmax = bwrite	!writes min/max in plot files
	berrintp = bverb	!writes errors in interpolation

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

	call iff_init_global_2d(nkn,nel,hkv,hev,date,time)

	call init_plot

	call shy_get_string_descriptions(id,nvar,ivars,strings)
	call choose_var(nvar,ivars,strings,varline,ivarplot,bvect)

	bsect = getisec() /= 0
	call setlev(layer)

	if( binfo ) stop

	!--------------------------------------------------------------
	! initialize volume
	!--------------------------------------------------------------

	shy_zeta = 0.
	call shy_make_volume

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

	call init_regular
	call info_regular(bregplot,nx,ny,dx,dy)
	if( bregplot .and. .not. bquiet ) then
	  write(6,*) 'regular grid plotting: ',nx,ny,dx,dy
	end if

        call initialize_color

	call qopen

!--------------------------------------------------------------
! loop on data (time loop)
!--------------------------------------------------------------

	dtime = 0.
	cv3 = 0.
	cv3all = 0.

	do

	 !--------------------------------------------------------------
	 ! read new data set
	 !--------------------------------------------------------------

	 call read_records(id,dtime,nvar,nndim,nlvdi,idims
     +				,cv3,cv3all,ierr)

         if(ierr.ne.0) exit

	 call dts_convert_to_atime(datetime_elab,dtime,atime)

	 irec = irec + 1
	 nread = nread + nvar

	 !--------------------------------------------------------------
	 ! see if we are in time window
	 !--------------------------------------------------------------

	 bcycle = .false.
	 if( elabtime_over_time(atime) ) exit
	 if( .not. elabtime_in_time(atime) ) bcycle = .true.
	 if( ifreq > 0 .and. mod(irec,ifreq) /= 0 ) bcycle = .true.

	 if( bcycle ) then
	   if( bverb ) call shy_write_time2(irec,atime,0)
	   cycle
	 end if

	 call ptime_set_dtime(dtime)

	 call shy_make_zeta(ftype)
	 !call shy_make_volume		!comment for constant volume

	 !--------------------------------------------------------------
	 ! find out what to plot
	 !--------------------------------------------------------------

	 ivars(:) = idims(4,:)
	 call get_vars_to_plot(nvar,ivars,ivar3,ivarplot,bvect,ivnum,ivs)

	 ivar = ivarplot(1)
	 iv = ivs(1)
	 if( iv == 0 ) then
	   write(6,*) 'no such variable in time record: ',ivar,iv
	   stop 'error stop: no such variable'
	 end if

	 !--------------------------------------------------------------
	 ! extract variables
	 !--------------------------------------------------------------

	 if( bhydro ) then
	   znv(1:nkn) = cv3all(1,1:nkn,1)
	   zenv = reshape(cv3all(1,:,2),(/3,nel/))
	 end if

	 do i=1,2

	 iv = ivs(i)
	 if( iv == 0 ) cycle		!this happens for scalar
	 n = idims(1,iv)
	 m = idims(2,iv)
	 lmax = idims(3,iv)
	 ivar = idims(4,iv)
	 nn = n * m
	 if( m /= 1 ) stop 'error stop: m/= 1'
	 if( n /= nkn .and. n /= nel ) then
	      write(6,*) 'n,nkn,nel: ',n,nkn,nel
	      stop 'error stop: n'
	 end if

	 cv3(:,:) = cv3all(:,:,iv)
	 if( b2d ) then
	    call shy_make_vert_aver(idims(:,iv),nndim,cv3,cv2)
	 else if( lmax == 1 ) then
	   cv2(:) = cv3(1,:)
	 else
	   call extlev(layer,nlvdi,n,il,cv3,cv2)
	 end if

	 if( bvect ) then
	   call directional_insert(bvect,bregdata,ivar,ivar3
     +					,ivarplot,n,cv2,ivel)
	    if( iv == 3 ) utlnv(:,1:nel) = cv3(:,1:nel)
	    if( iv == 4 ) vtlnv(:,1:nel) = cv3(:,1:nel)
	 end if

	 end do

	 call make_vertical_velocity
	 call extlev(layer,nlvdi+1,n,il,wauxv,wsnv)	!average over layer

	 !------------------------------------------
	 ! from here plotting
	 !------------------------------------------

	 if( .not. bvect ) ivel = 0
	 bvel = ivel > 0
	 if( bvect .and. .not. bvel ) stop 'error stop: ivel == 0'

	 nplot = nplot + 1

	 if( bsdebug ) then
	   write(6,*) 'plotting: ',ivar,layer,n,ivel
	   write(6,*) 'plotting: ',ivar3,ivarplot
	 end if
	 btime = bverb .or. ivar > 0
	 btime = btime .and. .not. bsilent
	 !write(6,*) btime,bverb,bsilent,ivar
	 if( btime ) then
	   call shy_write_time2(irec,atime,ivar)
	 end if

	 call make_mask(layer)

	 if( bsect ) then
	   call plot_sect(bvel,cv3)
	 else if( bvel ) then
	   call plovect(ivel,'3D ',bregdata)
	 else
           call plo_scal_val(n,cv2,varline)
	 end if

	 !------------------------------------------
	 ! end of plotting
	 !------------------------------------------

	end do		!time do loop

!--------------------------------------------------------------
! end of loop on data (time loop)
!--------------------------------------------------------------

	call qclose

	if( .not. bsilent ) then
	  write(6,*) 'total number of plots: ',nplot
	end if

!--------------------------------------------------------------
! final write of variables
!--------------------------------------------------------------

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   71	continue
	write(6,*) 'ftype = ',ftype,'  nvar = ',nvar
	write(6,*) 'nvar should be 4'
	stop 'error stop shyelab: ftype,nvar'
   74	continue
	stop 'error stop shyelab: general error...'
   75	continue
	write(6,*) 'error writing header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: writing header'
   76	continue
	write(6,*) 'ftype = ',ftype,'  expecting 1 or 2'
	stop 'error stop shyelab: ftype'
   77	continue
	write(6,*) 'error reading header, ierr = ',ierr
	write(6,*) 'file = ',trim(file)
	stop 'error stop shyelab: reading header'
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine plot_fem_file

	use clo
	!use elabutil
	use plotutil
	use elabtime
	use shyfile
	use shyutil

        use basin
        use levels
        use evgeom
        use mod_depth
        use mod_geom
        use mod_hydro_plot

	implicit none

	logical bhasbasin,bregdata,bvect,bvel
	logical bsect,bskip
	logical bintp
	logical bregplot
	integer nx,ny
	integer i,ierr,iformat,irec,l,nplot,ivel,ivar
	integer isphe,iunit,lmax,lmax0,np,np0,npaux
	integer ntype,nvar,nvar0,nvers
	integer date,time
	integer datetime(2)
	integer itype(2)
	integer ivarplot(2),ivs(2)
	real x0,y0,dx,dy
	real regpar(7)
	real flag
	double precision dtime,atime,atime0
	character*80 line,string,varline
        real,allocatable :: data2d(:)
        real,allocatable :: data3d(:,:)
        real,allocatable :: data2ddir(:)
        real,allocatable :: data3ddir(:,:)
        real,allocatable :: data(:,:,:)
        real,allocatable :: dext(:)
        real,allocatable :: hd(:)
        integer,allocatable :: il(:)
        !real,allocatable :: hlv(:)
        !integer,allocatable :: ilhkv(:)
        integer,allocatable :: ivars(:)
        character*80, allocatable :: strings(:)

	integer getisec
	real getpar

        !--------------------------------------------------------------
        ! set command line parameters
        !--------------------------------------------------------------

        call init_nls_fnm
        call read_str_files(-1)
        call read_str_files(ivar3)

        !--------------------------------------------------------------
        ! open input files
        !--------------------------------------------------------------

	infile = femfilename
        if( infile .eq. ' ' ) stop

        np = 0
        call fem_file_read_open(infile,np,iformat,iunit)
        if( iunit .le. 0 ) stop

        call fem_file_get_format_description(iformat,line)

        !--------------------------------------------------------------
        ! set up params and modules
        !--------------------------------------------------------------

	!--------------------------------------------------------------
	! read first record
	!--------------------------------------------------------------

        call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

        if( ierr .ne. 0 ) goto 99

        if( .not. bquiet ) then
          write(6,*) 'file name: ',infile(1:len_trim(infile))
          write(6,*) 'format: ',iformat,"  (",line(1:len_trim(line)),")"
          write(6,*) 'nvers:  ',nvers
          write(6,*) 'np:     ',np
          write(6,*) 'lmax:   ',lmax
          write(6,*) 'nvar:   ',nvar
          write(6,*) 'ntype:  ',ntype
        end if

	nlv = lmax
        call levels_init(np,np,nlv)	!first call - will be changed later
        call fem_file_make_type(ntype,2,itype)

        call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)
        if( ierr .ne. 0 ) goto 98

        if( lmax > 1 .and. .not. bquiet ) then
          write(6,*) 'vertical layers: ',lmax
          write(6,*) hlv
        end if
        if( itype(1) .gt. 0 .and. .not. bquiet ) then
          write(6,*) 'date and time: ',datetime
        end if
	bregdata = ( itype(2) > 0 )
	if( bregdata ) dflag = regpar(7)
        if( bregdata .and. .not. bquiet ) then
          write(6,*) 'regpar: '
	  call printreg(regpar)
        end if

        !--------------------------------------------------------------
        ! configure basin
        !--------------------------------------------------------------

	bhasbasin = basfilename /= ' '

	!write(6,*) 'bhasbasin,bregdata: ',bhasbasin,bregdata
	call mod_hydro_set_regpar(regpar)
	if( bhasbasin ) then
          call read_command_line_file(basfilename)
	else if( bregdata ) then
	  call bas_insert_regular(regpar)
	else	!should not be possible
	  write(6,*) 'the fem file is not a regular file'
	  write(6,*) 'the parameters are given on an unstructured grid'
	  write(6,*) 'in order to plot them basin file .bas is needed'
	  write(6,*) 'please specify the .bas file on the command line'
	  write(6,*) 'bhasbasin,bregdata: ',bhasbasin,bregdata
	  stop 'error stop plot_fem_file: need basin'
	end if

	if( bhasbasin .or. bregdata ) then
	  call bash_verbose(bsdebug)
	  call ev_set_verbose(.not.bquiet)
          call ev_init(nel)
          isphe = nint(getpar('isphe'))
          call set_coords_ev(isphe)
          call set_ev
          call get_coords_ev(isphe)
          call putpar('isphe',float(isphe))

          call mod_geom_init(nkn,nel,ngr)
          call set_geom

          call levels_init(nkn,nel,nlv)
          call mod_depth_init(nkn,nel)

	  npaux = nel
	  if( bregdata .and. np > nkn ) then
	    write(6,*) 'changing dims: ',npaux,np,nkn
	    call mod_hydro_plot_init(np,2*np,nlv,np)
	  else
	    call mod_hydro_plot_init(nkn,nel,nlv,npaux)
	  end if

          call makehev(hev)
          call makehkv(hkv)
          call allocate_2d_arrays
	end if

	if( bregdata ) then		!data is in regular format
	  bintp = .false.
	  if( bhasbasin ) bintp = .true.
	  if( bregall ) bintp = .false.
	else			!data is 1 value, BC or on nodes
	  bintp = .true.
	  if( bhasbasin ) then
	    if( nkn /= np ) goto 93
	    if( basintype /= 'bas' ) goto 92
	  else
	    !goto 94
	  end if
	end if

        nvar0 = nvar
        lmax0 = lmax
        np0 = np
        allocate(strings(nvar))
        allocate(ivars(nvar))
        allocate(dext(nvar))
        allocate(data2d(np))
        allocate(data3d(lmax,np))
        allocate(data2ddir(np))
        allocate(data3ddir(lmax,np))
        allocate(data(lmax,np,nvar))
        allocate(hd(np))
        allocate(il(np))

	!--------------------------------------------------------------
	! choose variable to plot
	!--------------------------------------------------------------

        do i=1,nvar
          call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
          if( ierr .ne. 0 ) goto 97
          strings(i) = string
	  call string2ivar(strings(i),ivars(i))
        end do

	call choose_var(nvar,ivars,strings,varline,ivarplot,bvect)
	call get_vars_to_plot(nvar,ivars,ivar3,ivarplot,bvect,ivnum,ivs)
	call set_ivel(ivar3,ivel)
	bvel = ivel > 0
	if( bvect .and. .not. bvel ) stop 'error stop: ivel == 0'
        if( ivs(1) == 0 ) then
          write(6,*) 'no such variable in file: ',ivar3,ivarplot,ivs
	  stop 'error stop plot_fem_file'
	else if( ivs(2) == 0 .and. bvect ) then
          write(6,*) 'need two variables for directional plot ',ivs
	  stop 'error stop plot_fem_file'
	end if
	if( bsdebug ) write(6,*) 'what to plot: ',ivar3,ivarplot,ivs

	if( binfo ) stop

	bminmax = bwrite	!writes min/max in plot files
	berrintp = bverb	!writes errors in interpolation

	if( .not. bregdata .and. .not. bhasbasin ) goto 94

	!--------------------------------------------------------------
	! close and re-open file
	!--------------------------------------------------------------

        close(iunit)

        np = 0
        call fem_file_read_open(infile,np,iformat,iunit)
        if( iunit .le. 0 ) stop

        !--------------------------------------------------------------
        ! set time
        !--------------------------------------------------------------

        call dts_convert_to_atime(datetime,dtime,atime)
        atime0 = atime          !absolute time of first record

	date = datetime(1)
	time = datetime(2)
        call ptime_init
        call ptime_set_date_time(date,time)
        call elabtime_date_and_time(date,time)
        call elabtime_set_minmax(stmin,stmax)
        call elabtime_set_inclusive(.false.)

	!--------------------------------------------------------------
	! initialize plot
	!--------------------------------------------------------------

        call init_plot

	bsect = getisec() /= 0
	call setlev(layer)
	b2d = layer == 0

	call init_regular
	call info_regular(bregplot,nx,ny,dx,dy)
	if( bregplot ) then
	  write(6,*) 'regular grid plotting: ',nx,ny,dx,dy
	end if

        call initialize_color

	call qopen

        !--------------------------------------------------------------
        ! loop on records
        !--------------------------------------------------------------

        irec = 0
	nplot = 0

        do
          !------------------------------------------------------------
          ! read headers of record and elab time
          !------------------------------------------------------------

          call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

          if( ierr .lt. 0 ) exit
          irec = irec + 1
          if( ierr .gt. 0 ) goto 99
          if( nvar .ne. nvar0 ) goto 96
          if( lmax .ne. lmax0 ) goto 96
          if( np .ne. np0 ) goto 96

          call dts_convert_to_atime(datetime,dtime,atime)
          call dts_format_abs_time(atime,line)
	  call ptime_set_atime(atime)

          if( bsdebug ) write(6,*) irec,atime,trim(line)

          call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)
          if( ierr .ne. 0 ) goto 98
	  call init_sigma_info(lmax,hlv)

	  bskip = .false.
	  if( elabtime_over_time(atime) ) exit
	  if( .not. elabtime_in_time(atime) ) bskip = .true.
	  if( ifreq > 0 .and. mod(irec,ifreq) /= 0 ) bskip = .true.

	  if( bskip ) then
	    if( bverb ) then
	      write(6,*) irec,atime,trim(line),'   ...skipping'
	    end if
	  else
	    if( .not. bsilent ) then
	      write(6,*) irec,atime,trim(line),'   ...plotting'
	    end if
	  end if

          !------------------------------------------------------------
          ! read data of record
          !------------------------------------------------------------

          do i=1,nvar
            if( bskip ) then
              call fem_file_skip_data(iformat,iunit
     +                          ,nvers,np,lmax,string,ierr)
            else
              call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,il,hd
     +                          ,lmax,data(1,1,i)
     +                          ,ierr)
            end if
            if( ierr .ne. 0 ) goto 97
            if( string .ne. strings(i) ) goto 95
          end do

	  if( bskip ) cycle

          !------------------------------------------------------------
          ! extract variable to plot
          !------------------------------------------------------------

	  nplot = nplot + 1
	  flag = dflag

	  data3d = data(:,:,ivs(1))
	  call adjust_levels_with_flag(nlvdi,np,il,flag,data3d)
	  if( bvect ) then
	    data3ddir = data(:,:,ivs(2))
	    call adjust_levels_with_flag(nlvdi,np,il,flag,data3ddir)
	  end if

          !------------------------------------------------------------
          ! extract or create layer to plot
          !------------------------------------------------------------

	  if( lmax == 1 ) then
	    data2d(:) = data3d(1,:)
	    if( bvect ) data2ddir(:) = data3ddir(1,:)
	  else if( b2d ) then
	    call fem_average_vertical(nlvdi,np,lmax,il,hlv,hd
     +					,data3d,data2d)
	    if( bvect ) then
	      call fem_average_vertical(nlvdi,np,lmax,il,hlv,hd
     +					,data3ddir,data2ddir)
	    end if
	  else
	    call extlev(layer,nlvdi,np,il,data3d,data2d)
	    if( bvect ) then
	      call extlev(layer,nlvdi,np,il,data3ddir,data2ddir)
	    end if
	  end if

	  !write(6,*) 'data: ',lmax,b2d,layer,(data2d(i),i=1,np,np/10)

          !------------------------------------------------------------
          ! plot variable(s)
          !------------------------------------------------------------

	  if( bvect ) then
	    ivar = ivars(ivs(1))
	    call directional_insert(bvect,bregdata,ivar,ivar3,ivarplot
     +					,np,data2d,ivel)
	    ivar = ivars(ivs(2))
	    call directional_insert(bvect,bregdata,ivar,ivar3,ivarplot
     +					,np,data2ddir,ivel)
	  end if

	  if( bregdata ) then
	    call reset_mask
	    if( bvect ) then
	      call plovect(ivel,'3D ',bregdata)
	    else
	      call ploreg(np,data2d,regpar,varline,bintp)
	    end if
	  else
            !call outfile_make_hkv(nkn,nel,nen3v,hm3v,hev,hkv)
	    if( np /= nkn ) stop 'cannot handle yet: internal error (9)'
	    ilhkv = il
            call ilhk2e(nkn,nel,nen3v,ilhkv,ilhv)
            call adjust_layer_index(nel,nlv,hev,hlv,ilhv)
	    call make_mask(layer)
	    if( bvect ) then
	      call plovect(ivel,'3D ',bregdata)
	    else
	      call ploval(np,data2d,varline)
	    end if
	  end if

	end do

        !--------------------------------------------------------------
        ! end loop on records
        !--------------------------------------------------------------

	call qclose

	if( .not. bsilent ) then
	  write(6,*) 'total number of plots: ',nplot
	end if

        !--------------------------------------------------------------
        ! end of routine
        !--------------------------------------------------------------

	return
   92   continue
        write(6,*) 'for non regular file we need bas, not grd file'
        stop 'error stop plot_fem_file: basin'
   93   continue
        write(6,*) 'incompatible node numbers: ',nkn,np
        stop 'error stop plot_fem_file: basin'
   94   continue
        write(6,*) 'fem file with non regular data needs basin'
        write(6,*) 'please specify basin on command line'
        stop 'error stop plot_fem_file: basin'
   95   continue
        write(6,*) 'strings not in same sequence: ',i
        write(6,*) string
        write(6,*) strings(i)
        stop 'error stop plot_fem_file: strings'
   96   continue
        write(6,*) 'nvar,nvar0: ',nvar,nvar0
        write(6,*) 'lmax,lmax0: ',lmax,lmax0    !this might be relaxed
        write(6,*) 'np,np0:     ',np,np0        !this might be relaxed
        write(6,*) 'cannot change number of variables'
        stop 'error stop plot_fem_file'
   97   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read data record of file'
        stop 'error stop plot_fem_file'
   98   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read second header of file'
        stop 'error stop plot_fem_file'
   99   continue
        write(6,*) 'record: ',irec
        write(6,*) 'cannot read header of file'
        stop 'error stop plot_fem_file'
	end

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine reset_mask

        call reset_dry_mask

	end

!***************************************************************

	subroutine make_mask(act_layer)

	use levels
        use mod_hydro_plot
        !use mod_hydro
        use plotutil

	implicit none

	integer act_layer

	integer level

	level = act_layer
	if( level == 0 ) level = 1

        call reset_dry_mask

        call set_level_mask(bwater,ilhv,level)        !element has this level
        call make_dry_node_mask(bwater,bkwater)       !copy elem to node mask

        call adjust_no_plot_area
        call make_dry_node_mask(bwater,bkwater)	      !copy elem to node mask
        if( bsdebug ) call info_dry_mask(bwater,bkwater)

	end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine initialize_color

	use plotutil

        implicit none

        integer icolor
        real getpar
	logical has_color_table

        call colsetup
	call admin_color_table

	if( has_color_table() ) call putpar('icolor',8.)

        icolor = nint(getpar('icolor'))
        call set_color_table( icolor )
        call set_default_color_table( icolor )

	if( bverb ) call write_color_table

        end

!***************************************************************

c*****************************************************************

        subroutine allocate_2d_arrays

        use mod_hydro_plot
        use mod_geom
        use mod_depth
        use evgeom
        use basin, only : nkn,nel,ngr,mbw
	use plotutil

        implicit none

        call ev_init(nel)
        call mod_geom_init(nkn,nel,ngr)

        call mod_depth_init(nkn,nel)

        if( bverb ) write(6,*) 'allocate_2d_arrays: ',nkn,nel,ngr

        end

c*****************************************************************

        subroutine read_command_line_file(file)

        use basin
        !use basutil
	use plotutil

        implicit none

        character*(*) file
        logical is_grd_file

        if( basin_is_basin(file) ) then
          if( .not. bsilent) write(6,*) 'reading BAS file ',trim(file)
          call basin_read(file)
          !breadbas = .true.
        else if( is_grd_file(file) ) then
          if( .not. bsilent) write(6,*) 'reading GRD file ',trim(file)
          call grd_read(file)
          call grd_to_basin
          call estimate_ngr(ngr)
	  call basin_set_read_basin(.true.)
          !breadbas = .false.
        else
          write(6,*) 'Cannot read this file: ',trim(file)
          stop 'error stop read_given_file: format not recognized'
        end if

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c routines for directional plots
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine set_ivel(ivar3,ivel)

! ivel == 3	wind
! ivel == 4	wave
! ivel == 5	velocities already prepared

	implicit none

	integer ivar3
	integer ivel

	if( ivar3 == 2 .or. ivar3 == 3 ) then
	  ivel = 5
	else if( ivar3 > 230 .and. ivar3 < 240 ) then
	  ivel = 4
	else if( ivar3 == 21 ) then
	  ivel = 3
	else
	  ivel = 0
	end if

	end

c*****************************************************************

	subroutine directional_init(nvar,ivars,ivar33,bvect,ivarplot)

! from ivar33 sets up bvect and ivarplot

	use plotutil

	implicit none

	integer nvar			!total number of variables in ivars
	integer ivars(nvar)		!id of variables in file
	integer ivar33			!id of what to plot
	logical bvect			!it is a vector variable (out)
	integer ivarplot(2)		!what variable id to use (out)

	logical bvel,bwave,bwind
	integer ivar,iv,nv,ivel

	ivarplot = ivar33

	bvel  = ivar33 == 2 .or. ivar33 == 3		!velocity/transport
	bwave = ivar33 > 230 .and. ivar33 < 240		!wave plot
	bwind = ivar33 == 21				!wind

	bvect = .false.
	bvect = bvect .or. bvel
	bvect = bvect .or. bwave
	bvect = bvect .or. bwind

	if( .not. bvect ) return

	if( bvel ) then
	  ivarplot = 2
	  if( any( ivars == 3 ) ) ivarplot = 3 	!transports in shy file
	end if
	if( bwave ) ivarplot = (/ivar33,233/)
	if( bwind ) ivarplot = 21

	nv = 0
	do iv=1,nvar
	  ivar = ivars(iv)
	  if( any( ivarplot == ivar ) ) then
	    nv = nv + 1
	  end if
	end do

	if( nv == 2 ) then
	  !write(6,*) 'can plot vector: ',ivarplot
	else if( nv == 1 ) then
	  if( bvect ) then
	    bvect = .false.
	    if( bverb ) write(6,*) 'can plot vector only as scalar: '
     +					,ivarplot
	  end if
	else if( nv == 0 ) then
	  write(6,*) '*** file does not contain needed varid: ',ivar33
          !stop 'error stop shyplot'
	else
	  write(6,*) '*** error in variables: ',ivar33,nv,nvar
	  write(6,*) ivars(:)
          stop 'error stop shyplot'
	end if

	end

c*****************************************************************

	subroutine directional_insert(bvect,bregdata,ivar,ivar3,ivarplot
     +					,n,cv2,ivel)

	use basin
	use mod_hydro_plot

	implicit none

	logical bvect		!can we plot a vector variable?
	logical bregdata	!did we read regular data?
	integer ivar		!variable id read
	integer ivar3		!variable id requested
	integer ivarplot(2)	!variable ids needed
	integer n		!total number of data points
	real cv2(n)		!data (2d)
	integer ivel		!on return indicates if and what to plot

! ivel == 1	velocities (from transports)
! ivel == 2	transports
! ivel == 3	wind
! ivel == 4	wave
! ivel == 5	velocities already prepared
!
! ivar == 2	velocities
! ivar == 3	transports

	logical bwave,bvel,bwind
	integer, save :: iarrow = 0

	ivel = 0
	bonelem = .false.
	bisreg = bregdata		!this is a global value
	bistrans = .false.

	if( .not. bvect ) return

	bvel  = ivar3 == 2 .or. ivar3 == 3		!velocity/transport
	bwave = ivar3 > 230 .and. ivar3 < 240		!wave plot
	bwind = ivar3 == 21				!wind

	if( bvel ) then
	  if( ivar == 3 ) then				!we read transports
	    if( n /= nel ) goto 99			!only for shy files
	    iarrow = iarrow + 1
	    bonelem = .true.
	    bistrans = .true.
	    if( iarrow == 1 ) utrans(1:nel) = cv2(1:nel)
	    if( iarrow == 2 ) vtrans(1:nel) = cv2(1:nel)
	  end if
	  if( ivar == 2 ) then				!we read velocities
	    iarrow = iarrow + 1
	    if( n == nel ) then
	      bonelem = .true.
	      if( iarrow == 1 ) uvelem(1:n) = cv2(1:n)
	      if( iarrow == 2 ) vvelem(1:n) = cv2(1:n)
	    else if( n == nkn ) then
	      if( iarrow == 1 ) uvnode(1:n) = cv2(1:n)
	      if( iarrow == 2 ) vvnode(1:n) = cv2(1:n)
	    else if( bregdata ) then
	      if( iarrow == 1 ) uvnode(1:n) = cv2(1:n)
	      if( iarrow == 2 ) vvnode(1:n) = cv2(1:n)
	    else
	      goto 99
	    end if
	  end if
	  if( iarrow == 2 ) then
	    ivel = ivar3 - 1
	    !call make_vertical_velocity
	    wsnv = 0.		!FIXME
	  end if
	end if

	if( bwave ) then
	  if( ivarplot(1) == ivar ) then
	    iarrow = iarrow + 1
	    uvspeed(1:n) = cv2(1:n)
	  else if( ivarplot(2) == ivar ) then
	    iarrow = iarrow + 1
	    uvdir(1:n) = cv2(1:n)
	  end if
	  if( iarrow == 2 ) then
	    ivel = 4
	    call polar2xy(n,uvspeed,uvdir,uvnode,vvnode)
	  end if
	end if

	if( bwind ) then
	  if( ivar == 21 ) then
	    iarrow = iarrow + 1
	    if( iarrow == 1 ) uvnode(1:n) = cv2(1:n)
	    if( iarrow == 2 ) vvnode(1:n) = cv2(1:n)
	  end if
	  if( iarrow == 2 ) then
	    ivel = 3
	    uvspeed = sqrt( uvnode**2 + vvnode**2 )
	  end if
	end if

	if( ivel > 0 ) iarrow = 0	!reset for next records

	!write(6,*) 'directional: ',bwind,bisreg

	return
   96	continue
	write(6,*) 'can read wind data only on nodes: ',n,nkn
	stop 'error stop directional_insert: wind data not on nodes'
   97	continue
	write(6,*) 'can read wave data only on nodes: ',n,nkn
	stop 'error stop directional_insert: wave data not on nodes'
   98	continue
	write(6,*) 'can read velocities only on nodes: ',n,nkn
	stop 'error stop directional_insert: velocities not on nodes'
   99	continue
	write(6,*) 'can read transports only on elements: ',n,nel
	stop 'error stop directional_insert: transports not on elems'
	end

c*****************************************************************

	subroutine choose_var(nvar,ivars,strings,varline,ivarplot,bvect)

	use plotutil
	use levels

! choses variable to be plotted
!
! either ivnum (sequential number) or ivar3 must have been set
! both of them are set on return
!
! cwarguments set on return:  varline, ivarplot
! global values set: ivar3, ivnum

	implicit none

	integer nvar
	integer ivars(nvar)
	character*80 strings(nvar)
	character*80 varline
	integer ivarplot(2)

	logical bvect
	integer nv,iv,ivar,isub
	character*80 string

!	---------------------------------------------------
!	write contents of file to terminal
!	---------------------------------------------------

	if( .not. bquiet ) then
	  call shy_print_descriptions(nvar,ivars,strings)
	end if

!	---------------------------------------------------
!	if sequential number is given, use this one
!	---------------------------------------------------

	if( ivnum > 0 ) then
	  if( ivnum > nvar ) then
	    write(6,*) 'ivnum too big for nvar'
	    write(6,*) 'ivnum,nvar: ',ivnum,nvar
            stop 'error stop shyplot'
	  end if
	  ivar3 = ivars(ivnum)
	end if

!	---------------------------------------------------
!	only one variable in file - plot this one
!	---------------------------------------------------

	if( nvar == 1 .and. ivar3 == 0 ) then
	  ivnum = 1
	  ivar3 = ivars(ivnum)
	end if
        if( ivar3 == 0 ) then
          write(6,*) '*** no variable given to be plotted: ',ivar3
	  write(6,*) 'Please provide one of the following: '
	  write(6,*) '  -varnum  -varid  -varname'
          stop 'error stop shyplot'
        end if

!	---------------------------------------------------
!	prepare for directional plot and read STR file
!	---------------------------------------------------

	call directional_init(nvar,ivars,ivar3,bvect,ivarplot)
	call read_str_files(ivar3)

!	---------------------------------------------------
!	see if ivar3 is in file
!	---------------------------------------------------

	call ivar2string(ivar3,string,isub)
	if( .not. bsilent ) then
	  write(6,*) 'variable to be plotted:  ',trim(string),ivar3
	end if
	nv = 0
	do iv=1,nvar
	  ivar = ivars(iv)
	  !write(6,'(2i10,4x,a)') iv,ivar,trim(strings(iv))
	  !if( ivar == ivar3 ) nv = nv + 1
	  if( ivar3 == ivar .or. any( ivarplot == ivar ) ) then
	    nv = nv + 1
	    if( ivnum == 0 ) ivnum = iv
	  end if
	end do

!	---------------------------------------------------
!	error check
!	---------------------------------------------------

	!if( nv == 0 .and. .not. bvect ) then
	if( nv == 0 ) then
	  write(6,*) 'no such variable in file: ',trim(string),ivar3
          stop 'error stop shyplot'
        end if

	if( layer > nlv ) then
          write(6,*) 'no such layer: ',layer
          write(6,*) 'maximum layer available: ',nlv
          stop 'error stop shyplot'
	end if

!	---------------------------------------------------
!	make varline and write to terminal
!	---------------------------------------------------

	call mkvarline(ivar3,varline)

	if( bsdebug ) then
	  write(6,*) 
	  write(6,*) 'debug information for plotting:'
	  write(6,*) 'varline: ',trim(varline)
	  write(6,*) 'ivnum: ',ivnum
	  write(6,*) 'ivar3: ',ivar3
	  write(6,*) 'layer: ',layer
	  write(6,*) 
	end if

!	---------------------------------------------------
!	end of routine
!	---------------------------------------------------

	end

c*****************************************************************

	subroutine fem_average_vertical(nlvddi,np,lmax,ilhkv,hlv,hd
     +					,data,data2d)

	implicit none

	integer nlvddi
	integer np,lmax
	integer ilhkv(np)
	real hlv(nlvddi)
	real hd(np)
	real data(nlvddi,np)
	real data2d(np)

	real hl(nlvddi)
	integer nsigma,k,lm,l,nlv
	real z,hsigma,h,hh,zeps
	double precision vacu,dacu

        call get_sigma_info(nlv,nsigma,hsigma)
	z = 0.
	zeps=0.01

	do k=1,np

	  lm = min(lmax,ilhkv(k))
	  if( lm <= 1 ) then
	    data2d(k) = data(1,k)
	    cycle
	  end if
	  h = hd(k)
	  if( h < -990. ) h = hlv(lm)
	  !if( h == -1. ) h = 1.
	  if( h+z<zeps ) z = zeps-h
          call get_layer_thickness(lm,nsigma,hsigma
     +                          ,z,h,hlv,hl)

	  vacu = 0.
	  dacu = 0.
	  do l=1,lm
	    hh = hl(l)
	    vacu = vacu + data(l,k)*hh
	    dacu = dacu + hh
	  end do
	  data2d(k) = vacu / dacu

	end do

	end

c*****************************************************************

	subroutine adjust_levels_with_flag(nlvddi,np,ilhkv,flag,data3d)

	implicit none

	integer nlvddi,np
	integer ilhkv(np)
	real flag
	real data3d(nlvddi,np)

	integer i

	do i=1,np
	  if( ilhkv(i) <= 0 ) then
	    ilhkv(i) = 1
	    data3d(1,i) = flag
	  end if
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine get_var_to_plot(nvar,ivars,ivar,iv)

	implicit none

	integer nvar
	integer ivars(nvar)
	integer ivar		!desired variable
	integer iv		!return

	do iv=1,nvar
	  if( ivars(iv) == ivar ) return
	end do

	iv = 0

	end

c*****************************************************************

	subroutine get_vars_to_plot(nvar,ivars,ivar,ivarplot
     +					,bvect,ivnum,ivs)

	implicit none

	integer nvar
	integer ivars(nvar)
	integer ivar		!desired variable
	integer ivarplot(2)	!desired variable for vector plot
	logical bvect
	integer ivnum
	integer ivs(2)		!return

	integer iv,ia

	ia = 1
	ivs = 0

	if( bvect ) then
	  do iv=1,nvar
	    if( ivars(iv) == ivarplot(ia) ) then
	      if( ia <= 2 ) ivs(ia) = iv
	      ia = ia + 1
	    end if
	  end do
	  if( ia /= 3 ) then
	    write(6,*) 'for vector plot need 2 variables to read: ',ia-1
	    ivs = 0
	  end if
	else
	  !do iv=1,nvar
	  !  if( ivars(iv) == ivar ) then
	  !    ivs(1) = iv
	  !    return
	  !  end if
	  !end do
	  ivs(1) = ivnum
	end if

	end

c*****************************************************************

        subroutine shy_write_time2(irec,atime,ivar)

        implicit none

	integer irec
        double precision atime
        integer ivar

        character*20 dline,extra

        dline = ' '
	extra = ' '
	if( ivar > 0 ) extra = '   ...plotting'

        call dts_format_abs_time(atime,dline)
	write(6,1000) irec,ivar,atime,'  ',trim(dline),extra
 1000	format(2i8,f16.2,a,a,a)

        end

c*****************************************************************

