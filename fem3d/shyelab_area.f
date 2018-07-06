
!***************************************************************
!
! these routines read line(s) defining area where averaging will be done
!
! old file format: xy
! new file format: grd, bnd, xy
!
! ieflag and ikflag are set and used in elab routines
!
!***************************************************************

        subroutine handle_area

! this works for more than one line

        use basin
        use elabutil

        implicit none

	integer is,ie,nn
        integer		  :: np
	real, allocatable :: xx(:),yy(:)
	integer, allocatable :: ifl(:)
	integer, allocatable :: ikf(:),ief(:)

        ieflag = 1
        ikflag = 1

	if( .not. barea ) return

	np = 0
        call read_all_lines(areafile,np,xx,yy,ifl)
        if( np <= 0 ) goto 99
        allocate(xx(np),yy(np),ifl(np))
        call read_all_lines(areafile,np,xx,yy,ifl)

        ieflag = -1
        ikflag = 0
        allocate(ikf(nkn),ief(nel))

	ie = 0
	do
	  is = ie + 1
	  call get_next_line(np,ifl,is,ie,nn)
	write(6,*) nn,is,ie
	  if( nn == 0 ) exit
          call check_elements(nn,xx(is:ie),yy(is:ie),ief,ikf)
	  where( ikf == 1 ) ikflag = 1
	  where( ief == 1 ) ieflag = 1
	  where( ief == 0 .and. ieflag < 0 ) ieflag = 0
	end do

	!call write_grd_with_flags('flags.grd',ikflag,ieflag)

	baverbas = .true.

	return
   99	continue
	write(6,*) 'error reading area file: ',trim(areafile)
	stop 'error stop handle_area: area file'
	end subroutine handle_area

!***************************************************************
!
! Read coordinates of points in line
! for n == 0 only checks how many nodes to read
! for n > 0 reads nodes into nodes() (error if n is too small)
!
!***************************************************************

	subroutine get_next_line(n,ifl,is,ie,nn)

	implicit none

	integer istart,n,is,ie,nn
	integer ifl(n)

	integer i

	nn = 0
	if( is > n ) return
	if( ifl(is) /= 1 ) goto 99
	ie = is

	do i=is+1,n
	  if( ifl(i) == 1 ) exit
	  ie = i
	end do

	nn = ie - is + 1

	return
   99	continue
	write(6,*) 'error in line file: ',is,ifl(is)
	stop 'error stop get_next_line: file structure'
	end

!***************************************************************

	subroutine write_grd_with_flags(file,ikf,ief)

	use basin

	implicit none

	character*(*) file
	integer ikf(nkn)
	integer ief(nel)

	integer k,ie,ii
	real x,y

	open(1,file=file,status='unknown',form='formatted')

	do k=1,nkn
	  x = xgv(k)
	  y = ygv(k)
	  write(1,1000) 1,k,ikf(k),x,y
	end do

	do ie=1,nel
	  x = xgv(k)
	  y = ygv(k)
	  write(1,2000) 2,ie,ief(ie),3,(nen3v(ii,ie),ii=1,3)
	end do

	close(1)

	return
 1000	format(i1,2i10,2f14.6)
 2000	format(i1,6i10)
	end

!***************************************************************

