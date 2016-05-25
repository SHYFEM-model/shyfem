c
c $Id: flxelab.f,v 1.8 2008-11-20 10:51:34 georg Exp $
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
c 05.10.2015    ggu     started flxelab
c
c**************************************************************

	subroutine flxelab

	use clo
	use elabutil
	use elabtime

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates nos file

	implicit none

	integer, allocatable :: iusplit(:)

	integer, allocatable :: kflux(:)
	integer, allocatable :: nlayers(:)
	real, allocatable :: fluxes(:,:,:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	integer nread,nelab,nrec,nin
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer ii,i,j,l,k,lmax,node
	integer ip,nb,naccum
	integer kfluxm
	integer idtflx,nlmax,nsect
	integer nscdi,nfxdi
	integer date,time
	character*80 title,name
	character*20 dline
	character*80 basnam,simnam
	real rnull
	real cmin,cmax,cmed,vtot
	real zmin,zmax,zmed
	real href,hzmin
	character*1 :: what(5) = (/'u','v','z','m','a'/)

	integer iapini,ifileo
	integer ifem_open_file

c--------------------------------------------------------------

	nread=0
	nelab=0
	nrec=0
	rnull=0.
	rnull=-1.
	bopen = .false.

c--------------------------------------------------------------
c set command line parameters
c--------------------------------------------------------------

	call elabutil_init('FLX')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	modeb = 2
	call ap_init(bask,modeb,0,0)

	call open_nos_type('.flx','old',nin)

        call flx_is_flx_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop flxelab: not a valid flx file'
        end if

	call flx_peek_header(nin,nvers,nsect,kfluxm,nlmax)
	nscdi = nsect
	nfxdi = kfluxm
	nlvdi = nlmax
	nlv = nlmax

	allocate(kflux(kfluxm))
	allocate(nlayers(nscdi))
	allocate(fluxes(0:nlvdi,3,nscdi))
	allocate(iusplit(nscdi))
	iusplit = 0


        call flx_read_header(nin,nscdi,nfxdi,nlvdi
     +                          ,nvers
     +                          ,nsect,kfluxm,idtflx,nlmax
     +                          ,kflux,nlayers
     +                          )

        write(6,*) 'nvers      : ',nvers
        write(6,*) 'nsect      : ',nsect
        write(6,*) 'kfluxm     : ',kfluxm
        write(6,*) 'idtflx     : ',idtflx
        write(6,*) 'nlmax      : ',nlmax
        write(6,*) 'kflux      : '
        write(6,*) (kflux(i),i=1,kfluxm)
        write(6,*) 'nlayers       : '
        write(6,*) (nlayers(i),i=1,nsect)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	!call nos_get_date(nin,date,time)
	date = 0
	time = 0
	call elabtime_date_and_time(date,time)

	!--------------------------------------------------------------
	! averaging
	!--------------------------------------------------------------

	!call elabutil_set_averaging(nvar)

	btrans = .false.
	if( btrans ) then
	  allocate(naccu(istep))
	  allocate(accum(nlvdi,nkn,istep))
	else
	  allocate(naccu(1))
	  allocate(accum(1,1,1))
	end if
	naccum = 0
	naccu = 0
	accum = 0.

	!write(6,*) 'mode: ',mode,ifreq,istep

	!--------------------------------------------------------------
	! open output file
	!--------------------------------------------------------------

	!iusplit = 0

	boutput = boutput .or. btrans
	bopen = boutput
	if( bsplit ) then
	  boutput = .false.
	  bopen = .false.
	end if

	if( bopen ) then
	  nb = ifileo(0,'out.flx','unform','new')
	  call flx_write_header(nb,nvers,nsect,kfluxm
     +			,idtflx,nlmax,kflux,nlayers)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it = 0
	if( .not. bquiet ) write(6,*)

	do

	 itold = it

	 call flx_read_record(nin,nvers,it,nlvdi,nsect,ivar
     +			,nlayers,fluxes,ierr)
         if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
         if(ierr.ne.0) exit
	 nread=nread+1

         if(ierr.ne.0) exit
	 nrec = nrec + 1

	 if( nrec == 1 ) itold = it
	 call flx_peek_record(nin,nvers,itnew,ierr)
	 !write(6,*) 'peek: ',it,itnew,ierr
	 if( ierr .ne. 0 ) itnew = it

	 if( .not. elabtime_check_time(it,itnew,itold) ) cycle

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( bwrite ) then
	    do l=1,nlv
	      do ii=1,3
	        !call mimar(xv(ii,:),kfluxm,zmin,zmax,rnull)
                !call aver(xv(ii,:),kfluxm,zmed,rnull)
	        !write(6,*) what(ii),' min,max,aver : ',zmin,zmax,zmed
	      end do
	      !do ii=1,kfluxm
	      !  s(ii) = sqrt( xv(1,ii)**2 + xv(2,ii)**2 )
	      !end do
	      !call mimar(s,kfluxm,zmin,zmax,rnull)
              !call aver(s,kfluxm,zmed,rnull)
	      !write(6,*) what(4),' min,max,aver : ',zmin,zmax,zmed
	    end do
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
	    !call make_vert_aver(nlvdi,nkn,ilhkv,cv3,vol3,cv2)
	  end if

	  if( bsplit ) then
            call split_flx(it,nlvdi,nsect,ivar,iusplit,fluxes)
            call fluxes_2d(it,nlvdi,nsect,ivar,fluxes)
            call fluxes_3d(it,nlvdi,nsect,ivar,nlayers,fluxes)
	  end if

	  if( boutput ) then
	    if( bverb ) write(6,*) 'writing to output: ',ivar
            call flx_write_record(nb,nvers,it
     +			,nlvdi,nsect,ivar,nlayers,fluxes)	!FIXME
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
	write(6,*) nread,' records read'
	!write(6,*) nrec ,' unique time records read'
	write(6,*) nelab,' records elaborated'
	write(6,*)

	if( bsplit ) then
	  write(6,*) 'output written to following files: '
	  write(6,*) '   p.[1-9]'
	else if( boutput ) then
	  write(6,*) 'output written to file out.flx'
	end if

	call ap_get_names(basnam,simnam)
	write(6,*) 'names used: '
	!write(6,*) 'basin: ',trim(basnam)
	write(6,*) 'simul: ',trim(simnam)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   85	continue
	write(6,*) 'it,itvar,i,ivar,nvar: ',it,itvar,i,ivar,nvar
	stop 'error stop noselab: time mismatch'
   92	continue
	write(6,*) 'incompatible basin: '
	write(6,*) 'nkn,nknnos: ',nkn,nknnos
	write(6,*) 'nel,nelnos: ',nel,nelnos
	stop 'error stop noselab: parameter mismatch'
   99	continue
	write(6,*) 'error writing to file unit: ',nb
	stop 'error stop noselab: write error'
	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine split_xv_0(iusplit,it,kfluxm,what,xv)

        implicit none

        integer iusplit(5,kfluxm)
	integer it
        integer kfluxm
        character*1 what(5)
        real xv(3,kfluxm)

	integer i,ii,iu
	real s
        character*80 name
        character*70 numb

	iu = 100

	if( iusplit(1,1) == 0 ) then
	  do i=1,kfluxm
	    do ii=1,5
	      iu = iu + 1
	      iusplit(ii,i) = iu
	      write(numb,'(i5)') i
	      numb = adjustl(numb)
	      name = what(ii) // '.' // numb
	      !write(6,*) 'opening file : ',iu,trim(name)
	      open(iu,file=name,form='formatted',status='unknown')
	    end do
	  end do
	end if

	do i=1,kfluxm
	  do ii=1,3
	    iu = iusplit(ii,i)
	    write(iu,*) it,xv(ii,i)
	  end do
	  iu = iusplit(4,i)
	  s = sqrt( xv(1,ii)**2 + xv(2,ii)**2 )
	  write(iu,*) it,s
	  iu = iusplit(5,i)
	  write(iu,*) it,xv(3,i),xv(1,i),xv(2,i),s
	end do

        end

c***************************************************************

        subroutine fluxes_2d(it,nlvddi,nsect,ivar,fluxes)

c writes 2d fluxes to file (only for ivar=0)

        implicit none

        integer it                      !time
        integer nlvddi                  !vertical dimension
        integer nsect                   !total number of sections
        integer ivar                    !type of variable (0: water fluxes)
        real fluxes(0:nlvddi,3,nsect)   !fluxes

        integer iunit,i,lmax,l,j
        real ptot(4,nsect)

        if( ivar .ne. 0 ) return

        do i=1,nsect
          ptot(1,i) = fluxes(0,1,i)                     !total
          ptot(2,i) = fluxes(0,2,i)                     !positive
          ptot(3,i) = fluxes(0,3,i)                     !negative
          ptot(4,i) = fluxes(0,2,i) + fluxes(0,3,i)     !absolute
        end do

        write(66,'(i10,20f10.2)') it,(ptot(1,i),i=1,nsect)      !total
        write(67,'(i10,20f10.2)') it,(ptot(2,i),i=1,nsect)      !positive
        write(68,'(i10,20f10.2)') it,(ptot(3,i),i=1,nsect)      !negative
        write(69,'(i10,20f10.2)') it,(ptot(4,i),i=1,nsect)      !absolute

c next is box format for Ali

        write(61,*) it
        write(61,*) 0
        write(61,*) nsect
        do i=1,nsect
          write(61,*) 0,0,(ptot(j,i),j=1,3)
        end do

        end

c****************************************************************

        subroutine fluxes_3d(it,nlvddi,nsect,ivar,nlayers,fluxes)

c writes 3d fluxes to file

        implicit none

        integer it                      !time
        integer nlvddi                  !vertical dimension
        integer nsect                   !total number of sections
        integer ivar                    !type of variable (0: water fluxes)
        integer nlayers(nsect)          !max layers for section
        real fluxes(0:nlvddi,3,nsect)   !fluxes

        integer iunit,i,lmax,l,j

        iunit = 160 + ivar

        write(iunit,*) it,nsect,ivar
        do i=1,nsect
          lmax = nlayers(i)
          write(iunit,*) i,lmax
          do l=0,lmax
            write(iunit,*) l,(fluxes(l,j,i),j=1,3)
          end do
        end do

        end

c****************************************************************

        subroutine split_flx(it,nlvddi,nsect,ivar,iusplit,fluxes)

c writes data to file name.number

        implicit none

        integer it
	integer nlvddi
	integer nsect
	integer ivar
        integer iusplit(nsect)
        real fluxes(0:nlvddi,3,nsect)   !fluxes

        integer i,ii,iu
        integer in
        character*80 numlin,file
        integer ialfa,ifileo

        if( ivar .ne. 0 ) return

	do i=1,nsect
	  iu = iusplit(i)
	  if( iu == 0 ) then
            in = ialfa(float(i),numlin,-1,-1)
            file = 'p' // '.' // numlin(1:in)
	    iu = ifileo(200,file,'form','new')
	    if( iu <= 0 ) then
	      write(6,*) 'cannot open file for writing: ',trim(file)
	      stop 'error stop split_flx: cannot open file'
	    end if
	    iusplit(i) = iu
	  end if

	  !write(iu,*) it,(fluxes(0,ii,i),ii=1,3)
	  write(iu,*) it,fluxes(0,1,i)
	end do

        end

c*******************************************************************


