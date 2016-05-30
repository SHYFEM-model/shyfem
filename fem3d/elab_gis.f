!
! utility routines for shyelab: elabutil
!
! revision log :
!
! 15.07.2015	ggu	written from scratch
! 22.09.2015	ggu	new routine open_shy_file()
! 10.10.2015	ggu	code added to handle FLX routines
! 22.02.2016	ggu	handle catmode
! 15.04.2016	ggu	handle gis files with substitution of colon
!
!************************************************************


c***************************************************************
c***************************************************************
c***************************************************************

        subroutine gis_write_record(dtime,ivar,np,nlvddi,il,cv,xv,yv)

c writes one record to file (3D)

        !use basin

        implicit none

	double precision dtime
        integer ivar,np,nlvddi
        integer il(np)
        real cv(nlvddi,np)
	real xv(np),yv(np)

	integer it
        integer i,l,lmax
	integer nout
        real x,y
	character*80 format,name
	character*20 line,dateline
	character*3 var

	integer ifileo

	it = nint(dtime)
	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)
	call i2s0(ivar,var)

	name = 'extract_'//var//'_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

        write(nout,*) it,np,ivar,dateline

	lmax = 1

        do i=1,np
          if( nlvddi > 1 ) lmax = il(i)
          x = xv(i)
          y = yv(i)
	  write(format,'(a,i5,a)') '(i10,2g14.6,i5,',lmax,'g14.6)'
          write(nout,format) i,x,y,lmax,(cv(l,i),l=1,lmax)
        end do

	close(nout)

        end

c***************************************************************

        subroutine gis_write_hydro(dtime,np,nlvddi,il,zv,uv,vv,xv,yv)

c writes one record to file (3D)

        !use basin

        implicit none

	double precision dtime
	integer np
        integer nlvddi
        integer il(np)
        real zv(np)
        real uv(nlvddi,np)
        real vv(nlvddi,np)
	real xv(np)
	real yv(np)

        integer l,lmax,nn,i,it
	integer nout
        real x,y
	character*80 format,name
	character*20 line,dateline
	character*3 var

	integer ifileo

	it = nint(dtime)
	call dtsgf(it,dateline)
	call gis_subst_colon(dateline,line)

	name = 'extract_hydro_'//line//'.gis'
        nout = ifileo(60,name,'form','new')
	!write(6,*) 'writing: ',trim(name)

        write(nout,*) it,np,0,dateline

	lmax = 1

        do i=1,np
          if( nlvddi > 1 ) lmax = il(i)
          x = xv(i)
          y = yv(i)
	  nn = 1 + 2*lmax
	  write(format,'(a,i5,a)') '(i10,2g14.6,i5,',nn,'g14.6)'
          write(nout,format) i,x,y,lmax,zv(i)
     +			,(uv(l,i),vv(l,i),l=1,lmax)
        end do

	close(nout)

        end

c***************************************************************

	subroutine gis_subst_colon(line_old,line_new)

	implicit none

	character*(*) line_old,line_new

	integer n,i

	n = min(len(line_old),len(line_new))
	line_new = line_old

	do i=1,n
	  if( line_new(i:i) == ':' ) line_new(i:i) = '_'
	  if( line_new(i:i) == ' ' ) line_new(i:i) = '_'
	end do

	end

c***************************************************************

        subroutine gis_write_connect

c writes connectivity

        use basin

        implicit none

	integer ie,ii

	open(1,file='connectivity.gis',form='formatted',status='unknown')

	write(1,*) nel
	do ie=1,nel
	  write(1,*) ie,(nen3v(ii,ie),ii=1,3)
	end do

	close(1)

	end

c***************************************************************
