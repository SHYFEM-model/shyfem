!
! $Id: subous.f,v 1.3 2004/09/28 13:31:55 georg Exp $
!
! utility routines to read/write OUS file - file type 161
!
! contents :
!
!        subroutine inious
!        subroutine setous(iunit,nvers,nkn,nel,nlv)
!        subroutine getous(iunit,nvers,nkn,nel,nlv)
!        subroutine delous(iunit)
!        subroutine dimous(iunit,nknddi,nelddi,nlvddi)
!
!        subroutine errous(iunit,routine,text)
!        subroutine findous_err(iunit,routine,text,n)
!        function findous(iunit)
!        subroutine infoous(iunit,iout)
!
!        subroutine ous_init(iunit,nversion)
!        subroutine ous_close(iunit)
!        subroutine ous_check_dimension(iunit,nknddi,nelddi,nlvddi)
!
!        subroutine ous_get_date(iunit,date,time)
!        subroutine ous_set_date(iunit,date,time)
!        subroutine ous_get_title(iunit,title)
!        subroutine ous_set_title(iunit,title)
!        subroutine ous_get_femver(iunit,femver)
!        subroutine ous_set_femver(iunit,femver)
!        subroutine ous_get_params(iunit,nkn,nel,nlv)
!        subroutine ous_set_params(iunit,nkn,nel,nlv)
!        subroutine ous_get_hparams(iunit,href,hzmin)
!        subroutine ous_set_hparams(iunit,href,hzmin)
!        subroutine ous_clone_params(iu_from,iu_to)
!
!	 subroutine ous_is_ous_file(iunit,nvers)
!
!        subroutine ous_read_header(iunit,nkn,nel,nlv,ierr)
!        subroutine ous_write_header(iunit,nkn,nel,nlv,ierr)
!        subroutine ous_read_header2(iu,ilhv,hlv,hev,ierr)
!        subroutine ous_write_header2(iunit,ilhv,hlv,hev,ierr)
!	 subroutine ous_read_record(iu,it,nlvddi,ilhv,z,ze,ut,vt,ierr)
!	 subroutine ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)
!
!        subroutine ous_next_record(iunit,it,ierr)
!        subroutine ous_back_record(iunit)
!        subroutine ous_skip_header(iunit,ierr)
!        subroutine ous_skip_record(iunit,it,ierr)
!
! revision log :
!
! 21.08.2003	ggu	version 1 implemented
! 01.09.2003	ggu	first revision
! 02.09.2003	ggu	last bug fixes (nvers=3 -> nvers=1)
! 22.09.2004	ggu	bug fix in rdous/wrous -> ie instead of k
! 08.06.2011	ggu	new routine delous(), check for end in read
! 18.01.2014	ggu	restructured, new date,time,femver
! 29.10.2014	ggu	new routine ous_is_ous_file()
!
! notes :
!
! Usage writing:
!
!       open file
!       call ous_init
!       call ous_set_title      (not obligatory)
!       call ous_set_date       (not obligatory)
!       call ous_set_femver     (not obligatory)
!       call ous_write_header
!       call ous_write_header2
!       call ous_write_record
!       ...
!       call ous_close
!
! Usage reading:
!
!       open file
!       call ous_init
!       call ous_read_header
!       call dimous
!       call ous_get_title      (not obligatory)
!       call ous_get_date       (not obligatory)
!       call ous_get_femver     (not obligatory)
!       call ous_read_header2
!       call ous_read_record
!       ...
!       call ous_close
!
! format of file:
!
! versions 1 and greater
!
!	ftype,nvers
!	nkn,nel,nlv
!	href,hzmin
!	title
!       date,time				(version 2)
!       femver					(version 2)
!
!	(ilhv(ie),ie=1,nel)
!	(hlv(l),l=1,nlv)
!	(hev(ie),ie=1,nel)
!	
!	it
!	(z(k),k=1,nkn)
!	(ze(i),i=1,3*nel)
!	((ut(l,ie),l=1,ilhv(k)),ie=1,nel)
!	((vt(l,ie),l=1,ilhv(k)),ie=1,nel)
!
!************************************************************
!************************************************************
!************************************************************
! internal routines
!************************************************************
!************************************************************
!************************************************************
!----------------------------------------------------------------------
        module ioous
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

	subroutine inious

! sets up initial common block - internal routine

	implicit none

	include 'ousinf.h'

	integer i,n

	logical binit
	save binit
	data binit /.false./

	if( binit ) return

	binit = .true.

	ousitem = 0

	do n=1,ndim
	  do i=0,nitdim
	    ousvar(i,n) = 0
	  end do
	  do i=1,nchdim
	    ouschar(i,n) = ' '
	  end do
	  do i=1,nrldim
	    ousreal(i,n) = 0.
	  end do
	end do

	end

!************************************************************

	subroutine setous(iunit,nvers,nkn,nel,nlv)

! sets up parameter common block - internal routine

	implicit none

	include 'ousinf.h'

	integer iunit,nvers,nkn,nel,nlv

	integer n

! we do not check if unit has already been opened -> open with ifileo
! changed -> before calling this ous_init has to be called

	n = findous(iunit)

        if( n .eq. 0 ) then	!yyyyyyyyyyyyyyyy
          n = findous(0)
        end if

        if( n .eq. 0 ) then
          call errous(iunit,'setous','Cannot find entry.')
        end if

	ousvar(0,n) = iunit	!yyyyyyyyyyyyyyyy
	if( nvers .gt. 0 ) ousvar(1,n) = nvers
	if(   nkn .gt. 0 ) ousvar(2,n) = nkn
	if(   nel .gt. 0 ) ousvar(3,n) = nel
	if(   nlv .gt. 0 ) ousvar(4,n) = nlv

	end

!************************************************************

	subroutine getous(iunit,nvers,nkn,nel,nlv)

! gets parameter common block - internal routine

	implicit none

	include 'ousinf.h'

	integer iunit,nvers,nkn,nel,nlv

	integer n

	n = findous(iunit)
        if( n .eq. 0 ) then
          call errous(iunit,'getous','File is not initialized.')
        end if

	nvers = ousvar(1,n)
	nkn   = ousvar(2,n)
	nel   = ousvar(3,n)
	nlv   = ousvar(4,n)

	end

!************************************************************

        subroutine delous(iunit)

! closes ous file internal structure - internal routine
!
! please note that the file has still to be closed manually

        implicit none

	include 'ousinf.h'

        integer iunit

        integer n,i

        call findous_err(iunit,'delous','File is not open, cannot close.',n)

        do i=0,nitdim
          ousvar(i,n) = 0
        end do
        do i=1,nchdim
          ouschar(i,n) = ' '
        end do
        do i=1,nrldim
          ousreal(i,n) = 0.
        end do

        end

!************************************************************

	subroutine dimous(iunit,nknddi,nelddi,nlvddi)

! checks dimension of arrays

	implicit none

	integer iunit,nknddi,nelddi,nlvddi

	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

        if( nkn .gt. nknddi ) goto 99
        if( nel .gt. nelddi ) goto 99
        if( nlv .gt. nlvddi ) goto 99

	return
   99   continue
        write(6,*) 'nkn,nknddi : ',nkn,nknddi
        write(6,*) 'nel,nelddi : ',nel,nelddi
        write(6,*) 'nlv,nlvddi : ',nlv,nlvddi
        stop 'error stop dimous: dimension error'
	end

!************************************************************
!************************************************************
!************************************************************

        subroutine errous(iunit,routine,text)

! error routine for ous - internal routine

        implicit none

        integer iunit
        character*(*) routine,text

        write(6,*) 'For unit ',iunit,' in routine ',routine
        write(6,*) text
        stop 'error stop errous'

        end

!************************************************************

        subroutine findous_err(iunit,routine,text,n)

! finds entry for iunit -> returns it in n or stops with error

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) routine,text
        integer n

        n = findous(iunit)

        if( n .eq. 0 ) then
          call errous(iunit,routine,text)
        end if

        end

!************************************************************

        function findous(iunit)

! finds entry for iunit - internal routine

        implicit none

        include 'ousinf.h'

        integer iunit
        integer findous

        integer n

        do n=1,min(ousitem+1,ndim)              !look at one entry more
          if( ousvar(0,n) .eq. iunit ) goto 1
        end do
        n = 0
    1   continue

        if( n .gt. ousitem ) ousitem = n
        findous = n

        end

!************************************************************

        subroutine infoous(iunit,iout)

! writes info for unit - internal routine

        implicit none

        include 'ousinf.h'

        integer iunit,iout

        integer n,i

        call findous_err(iunit,'ous_info','Cannot find entry.',n)

        write(iout,*) 'iunit = ',iunit,' position = ',n

	write(iout,*) 'integer'
        do i=0,nitdim
          write(iout,*) i,ousvar(i,n)
        end do
	write(iout,*) 'character'
        do i=1,nchdim
          write(iout,*) i,ouschar(i,n)
        end do
	write(iout,*) 'double precision'
        do i=1,nrldim
          write(iout,*) i,ousreal(i,n)
        end do

        end

!************************************************************
!************************************************************
!************************************************************
! public routines
!************************************************************
!************************************************************
!************************************************************

        subroutine ous_init(iunit,nversion)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nversion

        integer n,nvers

        call inious

        if( iunit .le. 0 ) then
          write(6,*) 'ous_init: Cannot initialize for this unit'
          write(6,*) 'iunit = ',iunit
          call errous(iunit,'ous_init','Impossible unit number.')
        end if

        nvers = nversion
        if( nvers .le. 0 ) nvers = maxvers

        if( nvers .gt. maxvers ) then
          write(6,*) 'ous_init: Impossible version number'
          write(6,*) 'nvers = ',nvers,'   maxvers = ',maxvers
          call errous(iunit,'ous_init','Impossible version number.')
        end if

        if( nvers .lt. maxcomp ) then
          write(6,*) 'ous_init: Old function call'
          write(6,*) 'nvers = ',nvers,'   maxcomp = ',maxcomp
          call errous(iunit,'ous_init','Old function call.')
        end if

	nvers = maxvers	!always write with highest version

        n = findous(iunit)
        if( n .ne. 0 ) then
          call errous(iunit,'ous_init','Unit already open.')
        end if

        n = findous(0)
        if( n .eq. 0 ) then
          call errous(iunit,'ous_init','No space left (ndim).')
        end if

        ousvar(0,n) = iunit
        ousvar(1,n) = nvers

        rewind(iunit)

        end

!************************************************************

        subroutine ous_close(iunit)

        implicit none

        integer iunit

        call delous(iunit)

        end

!************************************************************

        subroutine ous_check_dimension(iunit,nknddi,nelddi,nlvddi)

        implicit none

        integer iunit,nknddi,nelddi,nlvddi

        call dimous(iunit,nknddi,nelddi,nlvddi)

        end

!************************************************************
!************************************************************
!************************************************************

        subroutine ous_get_date(iunit,date,time)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer date,time

        integer n

        call findous_err(iunit,'ous_get_date','Cannot find entry.',n)

        date = ousvar(6,n)
        time = ousvar(7,n)

        end

!************************************************************

        subroutine ous_set_date(iunit,date,time)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer date,time

        integer n

        call findous_err(iunit,'ous_set_date','Cannot find entry.',n)

        ousvar(6,n) = date
        ousvar(7,n) = time

        end

!************************************************************

        subroutine ous_get_title(iunit,title)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) title

        integer n

        call findous_err(iunit,'ous_get_title','Cannot find entry.',n)

        title = ouschar(1,n)

        end

!************************************************************

        subroutine ous_set_title(iunit,title)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) title

        integer n

        call findous_err(iunit,'ous_set_title','Cannot find entry.',n)

        ouschar(1,n) = title

        end

!************************************************************

        subroutine ous_get_femver(iunit,femver)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) femver

        integer n

        call findous_err(iunit,'ous_get_femver','Cannot find entry.',n)

        femver = ouschar(2,n)

        end

!************************************************************

        subroutine ous_set_femver(iunit,femver)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) femver

        integer n

        call findous_err(iunit,'ous_set_femver','Cannot find entry.',n)

        ouschar(2,n) = femver

        end

!************************************************************

        subroutine ous_get_params(iunit,nkn,nel,nlv)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nkn,nel,nlv

        integer nvers

        call getous(iunit,nvers,nkn,nel,nlv)

        end

!************************************************************

        subroutine ous_set_params(iunit,nkn,nel,nlv)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nkn,nel,nlv

        call setous(iunit,0,nkn,nel,nlv)

        end

!************************************************************

        subroutine ous_get_hparams(iunit,href,hzmin)

        implicit none

        include 'ousinf.h'

        integer iunit
        double precision href,hzmin

        integer n

        call findous_err(iunit,'ous_get_hparams','Cannot find entry.',n)

	href  = ousreal(1,n)
	hzmin = ousreal(2,n)

	end

!************************************************************

        subroutine ous_set_hparams(iunit,href,hzmin)

        implicit none

        include 'ousinf.h'

        integer iunit
        double precision href,hzmin

        integer n

        call findous_err(iunit,'ous_set_hparams','Cannot find entry.',n)

	ousreal(1,n) = href
	ousreal(2,n) = hzmin

	end


!************************************************************

        subroutine ous_clone_params(iu_from,iu_to)

! clones data from one to other file
!
! second file must have already been opened and initialized with ous_init
! should be only used to write file -> nvers should be max version

        implicit none

        include 'ousinf.h'

        integer iu_from
        integer iu_to

        integer i,nf,nt

        call findous_err(iu_from,'ous_clone_params','Cannot find entry.',nf)
        call findous_err(iu_to,'ous_clone_params','Cannot find entry.',nt)

        do i=2,nitdim           !unit and version are not cloned
          ousvar(i,nt) = ousvar(i,nf)
        end do
        do i=1,nchdim
          ouschar(i,nt) = ouschar(i,nf)
        end do
        do i=1,nrldim
          ousreal(i,nt) = ousreal(i,nf)
        end do

        end

!************************************************************
!************************************************************
!************************************************************


        subroutine ous_is_ous_file(iunit,nvers)

! checks if iunit is open on ous file - returns nvers
!
! nvers == 0    no ous file (ntype is different) or read error
! nvers < 0     version number is wrong
! nvers > 0     good ous file

        implicit none

        include 'ousinf.h'

        integer iunit,nvers

        integer ntype

        nvers = 0
	if( iunit .le. 0 ) return

        read(iunit,end=1,err=1) ntype,nvers

        if( ntype .ne. ftype ) nvers = 0
        if( nvers .le. 0 .or. nvers .gt. maxvers ) nvers = -abs(nvers)

	rewind(iunit)

    1   continue

        end

!************************************************************
!************************************************************
!************************************************************

	subroutine ous_read_header(iunit,nkn,nel,nlv,ierr)

! before this ous_init has to be called

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

	integer iunit
	integer nkn,nel,nlv
	integer ierr

	integer n,nvers
	integer ntype,irec
	integer date,time
	double precision href,hzmin
	character*80 line
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

! initialize

	call inious

        call findous_err(iunit,'ous_read_header','Cannot find entry.',n)

! first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

! control version number and type of file

	if( ntype .ne. ftype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxvers ) goto 98

! next records

	date = 0
	time = 0

	irec = 2
	if( nvers .ge. 1 ) then
	  read(iunit,err=99)	 nkn,nel,nlv
	  read(iunit,err=99)	 href,hzmin
	  read(iunit,err=99)	 line
	else
	   stop 'error stop ous_read_header: internal error (1)'
	end if

	call setous(iunit,nvers,nkn,nel,nlv)
	call ous_set_hparams(iunit,href,hzmin)
	call ous_set_title(iunit,line)

	if( nvers .ge. 2 ) then
	  read(iunit,err=99)	 date,time
	  read(iunit,err=99)	 line
	  call ous_set_date(iunit,date,time)
	  call ous_set_femver(iunit,line)
	end if

	ierr=0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	return
   99	continue
	write(6,*) 'ous_read_header: Error encountered while'
	write(6,*) 'reading record number ',irec
	write(6,*) 'of OUS file header'
	write(6,*) 'nvers = ',nvers
	ierr=99
	return
   98	continue
	write(6,*) 'ous_read_header: Version not recognized : ',nvers
	ierr=98
	return
   97	continue
	write(6,*) 'ous_read_header: Wrong type of file : ',ntype
	write(6,*) 'Expected ',ftype
	ierr=97
	return
   91	continue
	write(6,*) 'ous_read_header: File is empty'
	ierr=91
	return
	end

!********************************************************************

	subroutine ous_write_header(iunit,nkn,nel,nlv,ierr)

! writes first record of OUS file

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

	integer iunit
	integer nkn,nel,nlv
	integer ierr

	integer n,nvers
	integer date,time
	double precision href,hzmin
	character*80 title,femver
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	call inious

	call findous_err(iunit,'ous_write_header','Cannot find entry.',n)

	nvers = maxvers
	call setous(iunit,nvers,nkn,nel,nlv)

        call ous_get_hparams(iunit,href,hzmin)
        call ous_get_title(iunit,title)
        call ous_get_date(iunit,date,time)
        call ous_get_femver(iunit,femver)

	write(iunit)		ftype,maxvers
	write(iunit)		nkn,nel,nlv
	write(iunit)		href,hzmin
	write(iunit)		title
        write(iunit)            date,time
        write(iunit)            femver

	ierr=0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	end

!************************************************************

	subroutine ous_read_header2(iu,ilhv,hlv,hev,ierr)

! reads second record of OUS file

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

	integer iu
	integer ilhv(*)
	double precision hlv(*)
	double precision hev(*)
        real,dimension(:),allocatable :: hl,he
	integer ierr

	logical bdata
	integer iunit
	integer l,ie,neli
	integer nvers,nkn,nel,nlv
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

        bdata = iu .gt. 0       !with negative unit number skip arrays
        iunit = abs(iu)

	call getous(iunit,nvers,nkn,nel,nlv)

#ifdef SINGLEP
        allocate(hl(nlv))
        allocate(he(nel))
#endif

	neli = nel

        if( .not. bdata ) then  !do not read arrays
	  neli = 0
          nel = 0
          nlv = 0
        else if( nlv .le. 1 ) then
          do ie=1,nel
            ilhv(ie) = 1
          end do
          hlv(1) = 10000.

          neli = 0
          nlv = 0
        end if

! read records

	if( nvers .ge. 1 ) then
	  read(iunit,err=99) (ilhv(ie),ie=1,neli)
#ifdef SINGLEP
	  read(iunit,err=99) (hl(l),l=1,nlv)
	  read(iunit,err=99) (he(ie),ie=1,nel)
          hlv(1:nlv)=dble(hl(1:nlv))
          hev(1:nkn)=dble(he(1:nkn))
#else
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
#endif
	else
	   stop 'error stop ous_read_header2: internal error (1)'
	end if

	ierr = 0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	return
   99	continue
	write(6,*) 'rsous: Error encountered while'
	write(6,*) 'reading second part of OUS file header'
	ierr=99
	return
	end

!************************************************************

	subroutine ous_write_header2(iunit,ilhv,hlv,hev,ierr)

! writes second record of OUS file

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

	integer iunit
	integer ilhv(*)
	double precision hlv(*)
	double precision hev(*)
        real,dimension(:),allocatable :: hl
        real,dimension(:),allocatable :: he
	integer ierr

	integer l,ie,neli
	integer nvers,nkn,nel,nlv
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	call getous(iunit,nvers,nkn,nel,nlv)

#ifdef SINGLEP
        allocate(hl(nlv))
        allocate(he(nel))
        hl=real(hlv(1:nlv))
        he=real(hev(1:nkn))
#endif

! only one layer

	neli = nel
        if( nlv .le. 1 ) then
          nlv = 0
          neli = 0
        end if

! write records

	write(iunit) (ilhv(ie),ie=1,neli)
#ifdef SINGLEP
	write(iunit) (hl(l),l=1,nlv)
	write(iunit) (he(ie),ie=1,nel)
#else
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)
#endif

	ierr = 0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	end

!************************************************************

	subroutine ous_read_record(iu,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

! reads data record of OUS file

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

	integer iu,it
	integer nlvddi
	integer ilhv(*)
	double precision z(*)
	double precision ze(*)
	double precision ut(nlvddi,*)
	double precision vt(nlvddi,*)
        real,dimension(:),allocatable :: rz,rze
        real,dimension(:,:),allocatable ::rut,rvt
	integer ierr

	integer l,k,ie,i,lmax
	integer nvers,nkn,nel,nlv
	integer iunit
	logical bdata
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

        bdata = iu .gt. 0       !with negative unit number only read header
        iunit = abs(iu)

	call getous(iunit,nvers,nkn,nel,nlv)

#ifdef SINGLEP
        allocate(rz(nkn))
        allocate(rze(3*nel))
        allocate(rut(nlv,nel))
        allocate(rvt(nlv,nel))
#endif

	if( bdata .and. nlvddi .lt. nlv ) goto 97

	if( .not. bdata ) then
	  nkn = 0
	  nel = 0
	end if

	if( nvers .ge. 1 ) then
	  read(iunit,end=88,err=98) it
#ifdef SINGLEP
	  read(iunit,end=99,err=99) (rz(k),k=1,nkn)
	  read(iunit,end=99,err=99) (rze(i),i=1,3*nel)
          z(1:nkn)=dble(rz(1:nkn))
          ze(1:3*nel)=dble(rze(1:3*nel))
#else
	  read(iunit,end=99,err=99) (z(k),k=1,nkn)
	  read(iunit,end=99,err=99) (ze(i),i=1,3*nel)
#endif
	  if( nlv .le. 1 ) then
#ifdef SINGLEP
	    read(iunit,end=99,err=99) (rut(1,ie),ie=1,nel)
	    read(iunit,end=99,err=99) (rvt(1,ie),ie=1,nel)
            ut(1,1:nel)=dble(rut(1,1:nel))
            vt(1,1:nel)=dble(rvt(1,1:nel))
#else
	    read(iunit,end=99,err=99) (ut(1,ie),ie=1,nel)
	    read(iunit,end=99,err=99) (vt(1,ie),ie=1,nel)
#endif
	  else
#ifdef SINGLEP
	    read(iunit,end=99,err=99) ((rut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	    read(iunit,end=99,err=99) ((rvt(l,ie),l=1,ilhv(ie)),ie=1,nel)
            ut(1:nlv,1:nel)=dble(rut(1:nlv,1:nel))
            vt(1:nlv,1:nel)=dble(rvt(1:nlv,1:nel))
#else
	    read(iunit,end=99,err=99) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	    read(iunit,end=99,err=99) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
#endif
	  end if
	else
	   stop 'error stop ous_read_record: internal error (1)'
	end if

	ierr=0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	return
   88	continue
	ierr=-1
	return
   97	continue
	write(6,*) 'ous_read_record:: nlvddi < nlv'
	write(6,*) 'nlvddi = ',nlvddi,'  nlv = ',nlv
	ierr=97
	return
   98	continue
	write(6,*) 'ous_read_record:: Error while reading'
	write(6,*) 'time record of OUS file'
	ierr=98
	return
   99	continue
	write(6,*) 'ous_read_record:: Error while reading'
	write(6,*) 'data record of OUS file'
	write(6,*) 'it = ',it
	ierr=99
	return
	end

!************************************************************

	subroutine ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

! writes data record of OUS file

	use shympi
        use timing

	implicit none

        include 'ousinf.h'

! arguments
	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	double precision z(*)
	double precision ze(*)
	double precision ut(nlvddi,*)
	double precision vt(nlvddi,*)
        real,dimension(:),allocatable :: rz
        real,dimension(:),allocatable :: rze
        real,dimension(:,:),allocatable :: rut
        real,dimension(:,:),allocatable :: rvt
	integer ierr
! local
	integer l,k,ie,i
	integer nvers,nkn,nel,nlv
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

	call getous(iunit,nvers,nkn,nel,nlv)
#ifdef SINGLEP
        allocate(rz(nkn))
        allocate(rze(3*nel))
        allocate(rut(nlvddi,nel))
        allocate(rvt(nlvddi,nel))
        rz=real(z(1:nkn))
        rze=real(ze(1:3*nel))
        rut=real(ut(1:nlvddi,1:nel))
        rvt=real(vt(1:nlvddi,1:nel))
#endif
	if( nlvddi .lt. nlv ) goto 97

	write(iunit) it
#ifdef SINGLEP
	write(iunit) (rz(k),k=1,nkn)
	write(iunit) (rze(i),i=1,3*nel)
#else
	write(iunit) (z(k),k=1,nkn)
	write(iunit) (ze(i),i=1,3*nel)
#endif
	if( nlv .le. 1 ) then
#ifdef SINGLEP
	  write(iunit) (rut(1,ie),ie=1,nel)
	  write(iunit) (rvt(1,ie),ie=1,nel)
#else
	  write(iunit) (ut(1,ie),ie=1,nel)
	  write(iunit) (vt(1,ie),ie=1,nel)
#endif
	else
#ifdef SINGLEP
	  write(iunit) ((rut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	  write(iunit) ((rvt(l,ie),l=1,ilhv(ie)),ie=1,nel)
#else
	  write(iunit) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	  write(iunit) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
#endif
	end if

	ierr=0

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

	return
   97	continue
	write(6,*) 'ous_write_record:: nlvddi < nlv'
	write(6,*) 'nlvddi = ',nlvddi,'  nlv = ',nlv
	ierr=97
	return
	end

!************************************************************

        subroutine ous_peek_record(iu,it,ierr)

! peeks into data record of OUS file

	use shympi
        use timing

        implicit none

! arguments
        integer iu,it
        integer ierr
! local
        integer l,k,lmax
        integer nvers,nkn,nel,nlv,nvar
        integer iunit,ios
        double precision time1

        if(ln_timing) time1 = shympi_wtime()

        iunit = abs(iu)

        call getous(iunit,nvers,nkn,nel,nlv)

        lmax = nlv

        if( nvers .ge. 1 ) then
           read(iunit,iostat=ios) it
        else
           write(6,*) 'nvers = ',nvers,'  iunit = ',iunit
           stop 'error stop ous_peek_record: internal error (1)'
        end if

        if( ios > 0 ) then
          write(6,*) 'ous_peek_record: Error while reading'
          write(6,*) 'time record of OUS file'
          ierr=98
          return
        end if

        backspace(iu)

        if( ios < 0 ) then
          ierr=-1
        else
          ierr=0
        end if

        if(ln_timing)  io_time = io_time + shympi_wtime() - time1

        end

!************************************************************
!************************************************************
!************************************************************

        subroutine ous_next_record(iunit,it,ierr)

! skips data record - only reads header of record

        implicit none

        integer iunit,it,ierr

        integer nlvddi
        integer ilhv(1)
	double precision z(1),ze(1)
        double precision ut(1,1),vt(1,1)

	nlvddi = 1
	call ous_read_record(-iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

        end

!************************************************************

        subroutine ous_back_record(iunit)

! skips back one data record (contains five reads)

        implicit none

        integer iunit

        backspace(iunit)
        backspace(iunit)
        backspace(iunit)
        backspace(iunit)
        backspace(iunit)

        end

!************************************************************

        subroutine ous_skip_header(iunit,ierr)

        implicit none

        integer iunit,ierr

        integer nkn,nel,nlv
        integer ilhv(1)
        double precision hlv(1)
        double precision hev(1)

        call ous_read_header(iunit,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) return
        call ous_read_header2(-iunit,ilhv,hlv,hev,ierr)

        end

!************************************************************

        subroutine ous_skip_record(iunit,it,ierr)

        implicit none

        integer iunit,it,ierr

        call ous_next_record(iunit,it,ierr)

        end

!************************************************************
!************************************************************
!************************************************************
! compatibility (old routine calls)
!************************************************************
!************************************************************
!************************************************************

	subroutine rfous(iunit,nvers,nkn,nel,nlv,href,hzmin,title,ierr)

	integer iunit,nvers
	integer nkn,nel,nlv
	double precision href,hzmin
	character*80 title
	integer ierr

	call ous_init(iunit,nvers)
	call ous_read_header(iunit,nkn,nel,nlv,ierr)
	call ous_get_title(iunit,title)
	call ous_get_hparams(iunit,href,hzmin)

	end

!************************************************************

	subroutine wfous(iunit,nvers,nkn,nel,nlv,href,hzmin,title,ierr)

	integer iunit,nvers
	integer nkn,nel,nlv
	double precision href,hzmin
	character*80 title
	integer ierr

	call ous_init(iunit,nvers)
	call ous_set_title(iunit,title)
	call ous_set_hparams(iunit,href,hzmin)
	call ous_write_header(iunit,nkn,nel,nlv,ierr)

	end

!************************************************************

	subroutine rsous(iunit,ilhv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhv(*)
	double precision hlv(*)
	double precision hev(*)
	integer ierr

	call ous_read_header2(iunit,ilhv,hlv,hev,ierr)

	end

!************************************************************

	subroutine wsous(iunit,ilhv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhv(*)
	double precision hlv(*)
	double precision hev(*)
	integer ierr

	call ous_write_header2(iunit,ilhv,hlv,hev,ierr)

	end

!************************************************************

	subroutine rdous(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	implicit none

	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	double precision z(*)
	double precision ze(*)
	double precision ut(nlvddi,*)
	double precision vt(nlvddi,*)
	integer ierr

	call ous_read_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	end

!************************************************************

	subroutine wrous(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

        use shympi
        use mpi_io_admin

	implicit none

	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	double precision z(*)
	double precision ze(*)
	double precision ut(nlvddi,*)
	double precision vt(nlvddi,*)
	integer ierr

        ierr=0
        if(shympi_partition_on_elements()) then
          call rebuild_structures(ilhv,z,ze,ut,vt)
          if(shympi_is_master()) then
             call ous_write_record(iunit,it,nlvddi,outIlhv,outZnv,outZenv,outUtlnv,outVtlnv,ierr)
          end if
        else
	  call ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)
        end if

	end

!************************************************************

!----------------------------------------------------------------------
        end module ioous
!----------------------------------------------------------------------
