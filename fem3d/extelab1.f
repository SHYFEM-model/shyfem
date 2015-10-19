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
c 05.10.2015    ggu     variables in xv were used in the wromg order - fixed
c
c**************************************************************

	subroutine extelab

	use clo
	use elabutil

        use basin
        use mod_depth
        use evgeom
        use levels

c elaborates nos file

	implicit none

	integer, allocatable :: iusplit(:,:)

	integer, allocatable :: knaus(:)
	real, allocatable :: hdep(:)
	real, allocatable :: xv(:,:)
	real, allocatable :: speed(:)

	integer, allocatable :: naccu(:)
	double precision, allocatable :: accum(:,:,:)

	integer nread,nelab,nrec,nin
	integer nvers
	integer nknnos,nelnos,nvar
	integer ierr
	integer it,ivar,itvar,itnew,itold,iaux
	integer ii,i,j,l,k,lmax,node
	integer ip,nb,naccum
	integer knausm
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

	call elabutil_init('EXT')

	!--------------------------------------------------------------
	! open input files
	!--------------------------------------------------------------

	modeb = 2
	call ap_init(bask,modeb,0,0)

	call open_nos_type('.ext','old',nin)

        call ext_is_ext_file(nin,nvers)
        if( nvers .le. 0 ) then
          write(6,*) 'nvers: ',nvers
          stop 'error stop extelab: not a valid ext file'
        end if

	call ext_peek_header(nin,nvers,knausm)
	nlv = 1

	allocate(knaus(knausm))
	allocate(hdep(knausm))
	allocate(xv(knausm,3))
	allocate(speed(knausm))
	allocate(iusplit(5,knausm))
	iusplit = 0

	call ext_read_header(nin,knausm,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,title)

        write(6,*) 'nvers      : ',nvers
        write(6,*) 'knausm     : ',knausm
        write(6,*) 'href,hzmin : ',href,hzmin
        write(6,*) 'title      : ',trim(title)
        write(6,*) 'knaus      : '
        write(6,*) (knaus(i),i=1,knausm)
        write(6,*) 'hdep       : '
        write(6,*) (hdep(i),i=1,knausm)

	!--------------------------------------------------------------
	! time management
	!--------------------------------------------------------------

	!call nos_get_date(nin,date,time)
	call elabutil_date_and_time

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
	  nb = ifileo(0,'out.ext','unform','new')
	  call ext_write_header(nb,knausm,nvers,knausm,knaus,hdep
     +                          ,href,hzmin,title)
	end if

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it = 0
	if( .not. bquiet ) write(6,*)

	do

	 itold = it

	 call ext_read_record(nin,nvers,it,knausm,xv,ierr)
         if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
         if(ierr.ne.0) exit
	 nread=nread+1

         if(ierr.ne.0) exit
	 nrec = nrec + 1

	 if( nrec == 1 ) itold = it
	 call ext_peek_record(nin,nvers,itnew,ierr)
	 !write(6,*) 'peek: ',it,itnew,ierr
	 if( ierr .ne. 0 ) itnew = it

	 if( .not. elabutil_check_time(it,itnew,itold) ) cycle

	  nelab=nelab+1

	  if( .not. bquiet ) then
	    dline = ' '
	    if( bdate ) call dtsgf(it,dline)
	    write(6,*) 'time : ',it,'  ',dline
	  end if

	  if( bwrite ) then
	    do l=1,nlv
	      do ii=1,3
	        call mimar(xv(:,ii),knausm,zmin,zmax,rnull)
                call aver(xv(:,ii),knausm,zmed,rnull)
	        write(6,*) what(ii),' min,max,aver : ',zmin,zmax,zmed
	      end do
	      do ii=1,knausm
	        speed(ii) = sqrt( xv(ii,1)**2 + xv(ii,2)**2 )
	      end do
	      call mimar(speed,knausm,zmin,zmax,rnull)
              call aver(speed,knausm,zmed,rnull)
	      write(6,*) what(4),' min,max,aver : ',zmin,zmax,zmed
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
            call split_xv(iusplit,it,knausm,what,xv)
	  end if

	  if( boutput ) then
	    if( bverb ) write(6,*) 'writing to output: ',ivar
            call ext_write_record(nb,nvers,it,knausm,xv)
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
	  write(6,*) '   [zuvma].[1-9]'
!	  do ivar=1,ndim
!	    if( iusplit(ivar) .gt. 0 ) then
!              write(name,'(i4)') ivar
!	      write(6,*) trim(adjustl(name))//'.nos'
!	    end if
!	  end do
	else if( boutput ) then
	  write(6,*) 'output written to file out.ext'
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

        subroutine split_xv(iusplit,it,knausm,what,xv)

        implicit none

        integer iusplit(5,knausm)
	integer it
        integer knausm
        character*1 what(5)
        real xv(knausm,3)

	integer i,ii,iu
	real speed
        character*80 name
        character*70 numb

	iu = 100

	if( iusplit(1,1) == 0 ) then
	  do i=1,knausm
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

	do i=1,knausm
	  do ii=1,3
	    iu = iusplit(ii,i)
	    write(iu,*) it,xv(i,ii)
	  end do
	  iu = iusplit(4,i)
	  speed = sqrt( xv(i,1)**2 + xv(i,2)**2 )
	  write(iu,*) it,speed
	  iu = iusplit(5,i)
	  write(iu,*) it,xv(i,3),xv(i,1),xv(i,2),speed
	end do

        end

c***************************************************************

