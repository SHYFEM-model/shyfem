c
c $Id: noselab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c 03.12.2001    ggu     some extra output -> place of min/max
c 09.12.2003    ggu     check for NaN introduced
c 07.03.2007    ggu     easier call
c 08.11.2008    ggu     do not compute min/max in non-existing layers
c 07.12.2010    ggu     write statistics on depth distribution (depth_stats)
c 06.05.2015    ggu     noselab started
c 05.06.2015    ggu     many more features added
c
c**************************************************************

	subroutine ouselab

	use clo
	use elabutil

	use basin
        use mod_depth
        use mod_hydro
        use mod_hydro_baro
        use evgeom
        use levels

c elaborates ous file

	implicit none

	real, allocatable :: u2v(:)
	real, allocatable :: v2v(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	real, allocatable :: hl(:)

	integer nread,nelab,nrec,nin
	integer nvers
	integer nknous,nelous
	integer invar,lmaxsp
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer i,l,k,lmax
	integer ip,nb,naccum
	character*80 title,name
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real zmin,zmax
	real umin,umax,vmin,vmax
	real volume

	integer iapini
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-999.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('OUS')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	call ap_init(bask,modeb,0,0)

        call open_ous_type('.ous','old',nin)

	call peek_ous_header(nin,nknous,nelous,nlv)

	if( bneedbasin ) then
	  if( nkn /= nknous .or. nel /= nelous ) goto 92
	else
	  nkn = nknous
	  nel = nelous
	end if

	call mod_hydro_baro_init(nel)
	call mod_depth_init(nkn,nel)
	call levels_init(nkn,nel,nlv)
	call mod_hydro_init(nkn,nel,nlv)

	allocate(u2v(nel),v2v(nel))
	allocate(hl(nlv))

	nlvdi = nlv
        call read_ous_header(nin,nkn,nel,nlvdi,ilhv,hlv,hev)
        call ous_get_params(nin,nkn,nel,nlv)

	call init_sigma_info(nlv,hlv)

	if( bneedbasin ) then
	  call outfile_make_hkv(nkn,nel,nen3v,hev,hkv)
	  call ev_init(nelous)
          call set_ev
	end if

	if( bverb ) call depth_stats(nkn,nlvdi,ilhkv)

	if( bnode ) then
	  call convert_internal_node(nodesp)
	  invar = 0
	  lmaxsp = ilhkv(nodesp)
	end if

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	call ous_get_date(nin,date,time)
	call elabutil_date_and_time

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	!call elabutil_set_averaging(nvar)
	btrans = .false.			!FIXME

	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	  naccum = 0
	  naccu = 0
	  accum = 0.
	end if

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	boutput = boutput .or. btrans
	bopen = boutput .and. .not. bsplit

	if( bopen ) then
          call open_ous_file('out','new',nb)
          call ous_init(nb,0)
          call ous_clone_params(nin,nb)
	  if( b2d ) then
	    call ous_set_params(nb,0,0,1)
	  end if
          call write_ous_header(nb,ilhv,hlv,hev)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	volume = 0.
	it = 0
	if( .not. bquiet ) write(6,*)

	do

	  itold = it

          call ous_read_record(nin,it,nlvdi,ilhv,znv,zenv
     +                          ,utlnv,vtlnv,ierr)

          if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
          if(ierr.ne.0) exit

	  nread=nread+1
	  nrec = nrec + 1

	  if( nrec == 1 ) itold = it
	  call ous_peek_record(nin,itnew,ierr)
	  if( ierr .ne. 0 ) itnew = it

	  if( .not. elabutil_check_time(it,itnew,itold) ) cycle

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( bwrite ) then

	    call mimar(znv,nkn,zmin,zmax,rnull)
            call comp_barotropic(nel,nlvdi,ilhv,utlnv,vtlnv,unv,vnv)
            call comp_vel2d(nel,hev,zenv,unv,vnv,u2v,v2v
     +                          ,umin,vmin,umax,vmax)
            !call compute_volume(nel,zenv,hev,volume)

            write(6,*) 'zmin/zmax : ',zmin,zmax
            write(6,*) 'umin/umax : ',umin,umax
            write(6,*) 'vmin/vmax : ',vmin,vmax
            write(6,*) 'volume    : ',volume

	  end if

	  if( btrans ) then
!	    call nos_time_aver(mode,i,ifreq,istep,nkn,nlvdi
!     +					,naccu,accum,cv3,boutput)
	  end if

	  if( baverbas ) then
!	    call make_aver(nlvdi,nkn,ilhkv,cv3,vol3
!     +                          ,cmin,cmax,cmed,vtot)
!	    call write_aver(it,ivar,cmin,cmax,cmed,vtot)
	  end if

	  if( b2d ) then
            call comp_barotropic(nel,nlvdi,ilhv,utlnv,vtlnv,unv,vnv)
	  end if

	  if( bsplit ) then
            !call get_split_iu(ndim,iusplit,ivar,nin,ilhkv,hlv,hev,nb)
	  end if

	  if( boutput ) then
	    if( bverb ) write(6,*) 'writing to output: ',ivar
	    if( bsumvar ) ivar = 30
	    if( b2d ) then
              call ous_write_record(nb,it,1,ilhv
     +			,znv,zenv,unv,vnv,ierr)
	    else
              call ous_write_record(nb,it,nlvdi,ilhv
     +			,znv,zenv,utlnv,vtlnv,ierr)
	    end if
            if( ierr .ne. 0 ) goto 99
	  end if

	  if( bnode ) then
	    ivar = 1
	    !call write_node(nodesp,nlvdi,cv3,it,ivar,lmaxsp)
	  end if

	end do		!time do loop

c--------------------------------------------------------------
c end of loop on data
c--------------------------------------------------------------

c--------------------------------------------------------------
c final write of variables
c--------------------------------------------------------------

	if( btrans ) then
	  !write(6,*) 'istep,naccu: ',istep,naccu
	  do ip=1,istep
	    naccum = naccu(ip)
	    !write(6,*) 'naccum: ',naccum
	    if( naccum > 0 ) then
	      write(6,*) 'final aver: ',ip,naccum
!	      call nos_time_aver(-mode,ip,ifreq,istep,nkn,nlvdi
!     +					,naccu,accum,cv3,boutput)
!	      if( bsumvar ) ivar = 30
!              call nos_write_record(nb,it,ivar,nlvdi,ilhkv,cv3,ierr)
              if( ierr .ne. 0 ) goto 99
	    end if
	  end do
	end if

c--------------------------------------------------------------
c write final message
c--------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' total records read'
	!write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
	write(6,*)

	if( boutput ) then
	  write(6,*) 'output written to file out.ous'
	end if

	call ap_get_names(basnam,simnam)
	write(6,*) 'names used: '
	write(6,*) 'basin: ',trim(basnam)
	write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknous: ',nkn,nknous
	write(6,*) 'nel,nelous: ',nel,nelous
	stop 'error stop ouselab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop ouselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

