
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c utility routines to read/write OUS file - file type 161
c
c contents :
c
c        subroutine inious
c        subroutine setous(iunit,nvers,nkn,nel,nlv)
c        subroutine getous(iunit,nvers,nkn,nel,nlv)
c        subroutine delous(iunit)
c        subroutine dimous(iunit,nknddi,nelddi,nlvddi)
c
c        subroutine errous(iunit,routine,text)
c        subroutine findous_err(iunit,routine,text,n)
c        function findous(iunit)
c        subroutine infoous(iunit,iout)
c
c        subroutine ous_init(iunit,nversion)
c        subroutine ous_close(iunit)
c        subroutine ous_check_dimension(iunit,nknddi,nelddi,nlvddi)
c
c        subroutine ous_get_date(iunit,date,time)
c        subroutine ous_set_date(iunit,date,time)
c        subroutine ous_get_title(iunit,title)
c        subroutine ous_set_title(iunit,title)
c        subroutine ous_get_femver(iunit,femver)
c        subroutine ous_set_femver(iunit,femver)
c        subroutine ous_get_params(iunit,nkn,nel,nlv)
c        subroutine ous_set_params(iunit,nkn,nel,nlv)
c        subroutine ous_get_hparams(iunit,href,hzmin)
c        subroutine ous_set_hparams(iunit,href,hzmin)
c        subroutine ous_clone_params(iu_from,iu_to)
c
c	 subroutine ous_is_ous_file(iunit,nvers)
c
c        subroutine ous_read_header(iunit,nkn,nel,nlv,ierr)
c        subroutine ous_write_header(iunit,nkn,nel,nlv,ierr)
c        subroutine ous_read_header2(iu,ilhv,hlv,hev,ierr)
c        subroutine ous_write_header2(iunit,ilhv,hlv,hev,ierr)
c	 subroutine ous_read_record(iu,it,nlvddi,ilhv,z,ze,ut,vt,ierr)
c	 subroutine ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)
c
c        subroutine ous_next_record(iunit,it,ierr)
c        subroutine ous_back_record(iunit)
c        subroutine ous_skip_header(iunit,ierr)
c        subroutine ous_skip_record(iunit,it,ierr)
c
c revision log :
c
c 21.08.2003	ggu	version 1 implemented
c 01.09.2003	ggu	first revision
c 02.09.2003	ggu	last bug fixes (nvers=3 -> nvers=1)
c 22.09.2004	ggu	bug fix in rdous/wrous -> ie instead of k
c 08.06.2011	ggu	new routine delous(), check for end in read
c 18.01.2014	ggu	restructured, new date,time,femver
c 29.10.2014	ggu	new routine ous_is_ous_file()
c
c notes :
c
c Usage writing:
c
c       open file
c       call ous_init
c       call ous_set_title      (not obligatory)
c       call ous_set_date       (not obligatory)
c       call ous_set_femver     (not obligatory)
c       call ous_write_header
c       call ous_write_header2
c       call ous_write_record
c       ...
c       call ous_close
c
c Usage reading:
c
c       open file
c       call ous_init
c       call ous_read_header
c       call dimous
c       call ous_get_title      (not obligatory)
c       call ous_get_date       (not obligatory)
c       call ous_get_femver     (not obligatory)
c       call ous_read_header2
c       call ous_read_record
c       ...
c       call ous_close
c
c format of file:
c
c versions 1 and greater
c
c	ftype,nvers
c	nkn,nel,nlv
c	href,hzmin
c	title
c       date,time				(version 2)
c       femver					(version 2)
c
c	(ilhv(ie),ie=1,nel)
c	(hlv(l),l=1,nlv)
c	(hev(ie),ie=1,nel)
c	
c	it
c	(z(k),k=1,nkn)
c	(ze(i),i=1,3*nel)
c	((ut(l,ie),l=1,ilhv(k)),ie=1,nel)
c	((vt(l,ie),l=1,ilhv(k)),ie=1,nel)
c
c************************************************************
c************************************************************
c************************************************************
c internal routines
c************************************************************
c************************************************************
c************************************************************

	subroutine inious

c sets up initial common block - internal routine

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

c************************************************************

	subroutine setous(iunit,nvers,nkn,nel,nlv)

c sets up parameter common block - internal routine

	implicit none

	include 'ousinf.h'

	integer iunit,nvers,nkn,nel,nlv

	integer n
	integer findous

c we do not check if unit has already been opened -> open with ifileo
c changed -> before calling this ous_init has to be called

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

c************************************************************

	subroutine getous(iunit,nvers,nkn,nel,nlv)

c gets parameter common block - internal routine

	implicit none

	include 'ousinf.h'

	integer iunit,nvers,nkn,nel,nlv

	integer n
	integer findous

	n = findous(iunit)
        if( n .eq. 0 ) then
          call errous(iunit,'getous','File is not initialized.')
        end if

	nvers = ousvar(1,n)
	nkn   = ousvar(2,n)
	nel   = ousvar(3,n)
	nlv   = ousvar(4,n)

	end

c************************************************************

        subroutine delous(iunit)

c closes ous file internal structure - internal routine
c
c please note that the file has still to be closed manually

        implicit none

	include 'ousinf.h'

        integer iunit

        integer n,i

        call findous_err(iunit,'delous'
     +                  ,'File is not open, cannot close.',n)

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

c************************************************************

	subroutine dimous(iunit,nknddi,nelddi,nlvddi)

c checks dimension of arrays

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

c************************************************************
c************************************************************
c************************************************************

        subroutine errous(iunit,routine,text)

c error routine for ous - internal routine

        implicit none

        integer iunit
        character*(*) routine,text

        write(6,*) 'For unit ',iunit,' in routine ',routine
        write(6,*) text
        stop 'error stop errous'

        end

c************************************************************

        subroutine findous_err(iunit,routine,text,n)

c finds entry for iunit -> returns it in n or stops with error

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) routine,text
        integer n

        integer findous

        n = findous(iunit)

        if( n .eq. 0 ) then
          call errous(iunit,routine,text)
        end if

        end

c************************************************************

        function findous(iunit)

c finds entry for iunit - internal routine

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

c************************************************************

        subroutine infoous(iunit,iout)

c writes info for unit - internal routine

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
	write(iout,*) 'real'
        do i=1,nrldim
          write(iout,*) i,ousreal(i,n)
        end do

        end

c************************************************************
c************************************************************
c************************************************************
c public routines
c************************************************************
c************************************************************
c************************************************************

        subroutine ous_init(iunit,nversion)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nversion

        integer n,nvers
        integer findous

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

c************************************************************

        subroutine ous_close(iunit)

        implicit none

        integer iunit

        call delous(iunit)

        end

c************************************************************

        subroutine ous_check_dimension(iunit,nknddi,nelddi,nlvddi)

        implicit none

        integer iunit,nknddi,nelddi,nlvddi

        call dimous(iunit,nknddi,nelddi,nlvddi)

        end

c************************************************************
c************************************************************
c************************************************************

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

c************************************************************

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

c************************************************************

        subroutine ous_get_title(iunit,title)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) title

        integer n

        call findous_err(iunit,'ous_get_title','Cannot find entry.',n)

        title = ouschar(1,n)

        end

c************************************************************

        subroutine ous_set_title(iunit,title)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) title

        integer n

        call findous_err(iunit,'ous_set_title','Cannot find entry.',n)

        ouschar(1,n) = title

        end

c************************************************************

        subroutine ous_get_femver(iunit,femver)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) femver

        integer n

        call findous_err(iunit,'ous_get_femver','Cannot find entry.',n)

        femver = ouschar(2,n)

        end

c************************************************************

        subroutine ous_set_femver(iunit,femver)

        implicit none

        include 'ousinf.h'

        integer iunit
        character*(*) femver

        integer n

        call findous_err(iunit,'ous_set_femver','Cannot find entry.',n)

        ouschar(2,n) = femver

        end

c************************************************************

        subroutine ous_get_params(iunit,nkn,nel,nlv)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nkn,nel,nlv

        integer nvers

        call getous(iunit,nvers,nkn,nel,nlv)

        end

c************************************************************

        subroutine ous_set_params(iunit,nkn,nel,nlv)

        implicit none

        include 'ousinf.h'

        integer iunit
        integer nkn,nel,nlv

        call setous(iunit,0,nkn,nel,nlv)

        end

c************************************************************

        subroutine ous_get_hparams(iunit,href,hzmin)

        implicit none

        include 'ousinf.h'

        integer iunit
        real href,hzmin

        integer n

        call findous_err(iunit,'ous_get_hparams','Cannot find entry.',n)

	href  = ousreal(1,n)
	hzmin = ousreal(2,n)

	end

c************************************************************

        subroutine ous_set_hparams(iunit,href,hzmin)

        implicit none

        include 'ousinf.h'

        integer iunit
        real href,hzmin

        integer n

        call findous_err(iunit,'ous_set_hparams','Cannot find entry.',n)

	ousreal(1,n) = href
	ousreal(2,n) = hzmin

	end


c************************************************************

        subroutine ous_clone_params(iu_from,iu_to)

c clones data from one to other file
c
c second file must have already been opened and initialized with ous_init
c should be only used to write file -> nvers should be max version

        implicit none

        include 'ousinf.h'

        integer iu_from
        integer iu_to

        integer i,nf,nt

        call findous_err(iu_from,'ous_clone_params'
     +                          ,'Cannot find entry.',nf)
        call findous_err(iu_to,'ous_clone_params'
     +                          ,'Cannot find entry.',nt)

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

c************************************************************
c************************************************************
c************************************************************


        subroutine ous_is_ous_file(iunit,nvers)

c checks if iunit is open on ous file - returns nvers
c
c nvers == 0    no ous file (ntype is different) or read error
c nvers < 0     version number is wrong
c nvers > 0     good ous file

        implicit none

        include 'ousinf.h'

        integer iunit,nvers

        integer ntype,ios

        nvers = 0
	if( iunit .le. 0 ) return

        read(iunit,iostat=ios) ntype,nvers
        if( ios /= 0 ) then
          nvers = 0
          return
        end if

        if( ntype .ne. ftype ) nvers = 0
        if( nvers .le. 0 .or. nvers .gt. maxvers ) nvers = -abs(nvers)

	rewind(iunit)

        end

c************************************************************
c************************************************************
c************************************************************

	subroutine ous_read_header(iunit,nkn,nel,nlv,ierr)

c before this ous_init has to be called

	implicit none

        include 'ousinf.h'

	integer iunit
	integer nkn,nel,nlv
	integer ierr

	integer n,nvers
	integer ntype,irec
	integer date,time
	real href,hzmin
	character*80 line

c initialize

	call inious

        call findous_err(iunit,'ous_read_header','Cannot find entry.',n)

c first record - find out what version

	irec = 1
	read(iunit,end=91,err=99) ntype,nvers

c control version number and type of file

	if( ntype .ne. ftype ) goto 97
	if( nvers .le. 0 .or. nvers .gt. maxvers ) goto 98

c next records

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
	backspace(iunit)
	ierr=91
	return
	end

c********************************************************************

	subroutine ous_write_header(iunit,nkn,nel,nlv,ierr)

c writes first record of OUS file

	implicit none

        include 'ousinf.h'

	integer iunit
	integer nkn,nel,nlv
	integer ierr

	integer n,nvers
	integer date,time
	real href,hzmin
	character*80 title,femver

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

	end

c************************************************************

	subroutine ous_read_header2(iu,ilhv,hlv,hev,ierr)

c reads second record of OUS file

	implicit none

        include 'ousinf.h'

	integer iu
	integer ilhv(*)
	real hlv(*)
	real hev(*)
	integer ierr

	logical bdata
	integer iunit
	integer l,ie,neli
	integer nvers,nkn,nel,nlv

        bdata = iu .gt. 0       !with negative unit number skip arrays
        iunit = abs(iu)

	call getous(iunit,nvers,nkn,nel,nlv)
	neli = nel

        if( .not .bdata ) then  !do not read arrays
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

c read records

	if( nvers .ge. 1 ) then
	  read(iunit,err=99) (ilhv(ie),ie=1,neli)
	  read(iunit,err=99) (hlv(l),l=1,nlv)
	  read(iunit,err=99) (hev(ie),ie=1,nel)
	else
	   stop 'error stop ous_read_header2: internal error (1)'
	end if

	ierr = 0

	return
   99	continue
	write(6,*) 'rsous: Error encountered while'
	write(6,*) 'reading second part of OUS file header'
	ierr=99
	return
	end

c************************************************************

	subroutine ous_write_header2(iunit,ilhv,hlv,hev,ierr)

c writes second record of OUS file

	implicit none

        include 'ousinf.h'

	integer iunit
	integer ilhv(*)
	real hlv(*)
	real hev(*)
	integer ierr

	integer l,ie,neli
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

c only one layer

	neli = nel
        if( nlv .le. 1 ) then
          nlv = 0
          neli = 0
        end if

c write records

	write(iunit) (ilhv(ie),ie=1,neli)
	write(iunit) (hlv(l),l=1,nlv)
	write(iunit) (hev(ie),ie=1,nel)

	ierr = 0

	end

c************************************************************

	subroutine ous_read_record(iu,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

c reads data record of OUS file

	implicit none

        include 'ousinf.h'

	integer iu,it
	integer nlvddi
	integer ilhv(*)
	real z(*)
	real ze(*)
	real ut(nlvddi,*)
	real vt(nlvddi,*)
	integer ierr

	integer l,k,ie,i,lmax
	integer nvers,nkn,nel,nlv
	integer iunit
	logical bdata

        bdata = iu .gt. 0       !with negative unit number only read header
        iunit = abs(iu)

	call getous(iunit,nvers,nkn,nel,nlv)

	if( bdata .and. nlvddi .lt. nlv ) goto 97

	if( .not. bdata ) then
	  nkn = 0
	  nel = 0
	end if

	if( nvers .ge. 1 ) then
	  read(iunit,end=88,err=98) it
	  read(iunit,end=99,err=99) (z(k),k=1,nkn)
	  read(iunit,end=99,err=99) (ze(i),i=1,3*nel)
	  if( nlv .le. 1 ) then
	    read(iunit,end=99,err=99) (ut(1,ie),ie=1,nel)
	    read(iunit,end=99,err=99) (vt(1,ie),ie=1,nel)
	  else
	    call linear2read(iunit,nlvddi,nel,ilhv,ut,ierr)
	    if( ierr /= 0 ) goto 99
	    call linear2read(iunit,nlvddi,nel,ilhv,vt,ierr)
	    if( ierr /= 0 ) goto 99
	    !read(iunit,end=99,err=99) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	    !read(iunit,end=99,err=99) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
	  end if
	else
	   stop 'error stop ous_read_record: internal error (1)'
	end if

	ierr=0

	return
   88	continue
	backspace(iunit)
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

c************************************************************

	subroutine ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

c writes data record of OUS file

	implicit none

        include 'ousinf.h'

c arguments
	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	real z(*)
	real ze(*)
	real ut(nlvddi,*)
	real vt(nlvddi,*)
	integer ierr
c local
	integer l,k,ie,i
	integer nvers,nkn,nel,nlv

	call getous(iunit,nvers,nkn,nel,nlv)

	if( nlvddi .lt. nlv ) goto 97

	write(iunit) it
	write(iunit) (z(k),k=1,nkn)
	write(iunit) (ze(i),i=1,3*nel)
	if( nlv .le. 1 ) then
	  write(iunit) (ut(1,ie),ie=1,nel)
	  write(iunit) (vt(1,ie),ie=1,nel)
	else
	  call linear2write(iunit,nlvddi,nel,ilhv,ut,ierr)
	  if( ierr /= 0 ) stop 'error stop ous_write_record'
	  call linear2write(iunit,nlvddi,nel,ilhv,vt,ierr)
	  if( ierr /= 0 ) stop 'error stop ous_write_record'
	  !write(iunit) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	  !write(iunit) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
	end if

	ierr=0

	return
   97	continue
	write(6,*) 'ous_write_record:: nlvddi < nlv'
	write(6,*) 'nlvddi = ',nlvddi,'  nlv = ',nlv
	ierr=97
	return
	end

c************************************************************

        subroutine ous_peek_record(iu,it,ierr)

c peeks into data record of OUS file

        implicit none

c arguments
        integer iu,it
        integer ierr
c local
        integer l,k,lmax
        integer nvers,nkn,nel,nlv,nvar
        integer iunit,ios

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

        end

c************************************************************
c************************************************************
c************************************************************

        subroutine ous_next_record(iunit,it,ierr)

c skips data record - only reads header of record

        implicit none

        integer iunit,it,ierr

        integer nlvddi
        integer ilhv(1)
	real z(1),ze(1)
        real ut(1,1),vt(1,1)

	nlvddi = 1
	call ous_read_record(-iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

        end

c************************************************************

        subroutine ous_back_record(iunit)

c skips back one data record (contains five reads)

        implicit none

        integer iunit

        backspace(iunit)
        backspace(iunit)
        backspace(iunit)
        backspace(iunit)
        backspace(iunit)

        end

c************************************************************

        subroutine ous_skip_header(iunit,ierr)

        implicit none

        integer iunit,ierr

        integer nkn,nel,nlv
        integer ilhv(1)
        real hlv(1)
        real hev(1)

        call ous_read_header(iunit,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) return
        call ous_read_header2(-iunit,ilhv,hlv,hev,ierr)

        end

c************************************************************

        subroutine ous_skip_record(iunit,it,ierr)

        implicit none

        integer iunit,it,ierr

        call ous_next_record(iunit,it,ierr)

        end

c************************************************************
c************************************************************
c************************************************************
c compatibility (old routine calls)
c************************************************************
c************************************************************
c************************************************************

	subroutine rfous	(iunit,nvers
     +				,nkn,nel,nlv
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

	integer iunit,nvers
	integer nkn,nel,nlv
	real href,hzmin
	character*80 title
	integer ierr

	call ous_init(iunit,nvers)
	call ous_read_header(iunit,nkn,nel,nlv,ierr)
	call ous_get_title(iunit,title)
	call ous_get_hparams(iunit,href,hzmin)

	end

c************************************************************

	subroutine wfous	(iunit,nvers
     +				,nkn,nel,nlv
     +				,href,hzmin
     +				,title
     +				,ierr
     +				)

	integer iunit,nvers
	integer nkn,nel,nlv
	real href,hzmin
	character*80 title
	integer ierr

	call ous_init(iunit,nvers)
	call ous_set_title(iunit,title)
	call ous_set_hparams(iunit,href,hzmin)
	call ous_write_header(iunit,nkn,nel,nlv,ierr)

	end

c************************************************************

	subroutine rsous(iunit,ilhv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhv(*)
	real hlv(*)
	real hev(*)
	integer ierr

	call ous_read_header2(iunit,ilhv,hlv,hev,ierr)

	end

c************************************************************

	subroutine wsous(iunit,ilhv,hlv,hev,ierr)

	implicit none

	integer iunit
	integer ilhv(*)
	real hlv(*)
	real hev(*)
	integer ierr

	call ous_write_header2(iunit,ilhv,hlv,hev,ierr)

	end

c************************************************************

	subroutine rdous(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	implicit none

	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	real z(*)
	real ze(*)
	real ut(nlvddi,*)
	real vt(nlvddi,*)
	integer ierr

	call ous_read_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	end

c************************************************************

	subroutine wrous(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	implicit none

	integer iunit,it
	integer nlvddi
	integer ilhv(*)
	real z(*)
	real ze(*)
	real ut(nlvddi,*)
	real vt(nlvddi,*)
	integer ierr

	call ous_write_record(iunit,it,nlvddi,ilhv,z,ze,ut,vt,ierr)

	end

c************************************************************

