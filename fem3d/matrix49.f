c
c reads/writes unit 49
c
c*******************************************************************
c
c Format of 2D files (with indices)
c
c Header:
c
c        read(49,*) nx,ny,nmax
c        read(49,*) (iboxes(1,k),iboxes(2,k),k=1,nmax)
c        read(49,*) (hmed(k),k=1,nmax)
c        read(49,*) (area(k),k=1,nmax)
c        read(49,*) (dxside(k),k=1,nmax)         
c        read(49,*) (dyside(k),k=1,nmax)
c        read(49,*) (sidex(k),k=1,nmax)         
c        read(49,*) (sidey(k),k=1,nmax)
c
c Record:
c
c        read(49,*,end=2) it
c        read(49,*) (zeta(k),k=1,nmax)
c        read(49,*) (uflux(k),k=1,nmax)
c        read(49,*) (vflux(k),k=1,nmax)
c
c*******************************************************************
c
c nx,ny           dimension of matrices
c nmax            dimension of pointer field
c 
c iboxes          pointer from array into matrix (only for debug purposes)
c 
c hmed            average depth of boxes
c area            area of boxes
c dxside	  equivalent dx of box (dxside,dyside are defined per box)
c dyside	  equivalent dy of box (dxside*dyside = area)
c sidex           length of vertical sides of boxes
c sidey           length of horizontal sides of boxes
c                 (format is the same as for fluxes afterwards)
c 
c it              time in seconds
c zeta            water level of boxes
c uflux           flux in x-direction
c vflux           flux in y-direction
c
c*******************************************************************

	subroutine wf49(nx,ny)

c writes header of unit 49

	implicit none

	integer nx,ny

	include 'regflux.h'

	integer k

	write(49,*) nx,ny,indhmx

	write(49,*) (iboxes(1,k),iboxes(2,k),k=1,indhmx)

	call mx2indx(1,1,nxdim,nydim,nx,ny,indhmx,hmed,value)
	call mx2indx(1,1,nxdim,nydim,nx,ny,indhmx,areav,value)
	call mx2indx(1,1,nxdim,nydim,nx,ny,indhmx,dxside,value)
	call mx2indx(1,1,nxdim,nydim,nx,ny,indhmx,dyside,value)
	call mx2indx(0,0,nxdim,nydim,nx,ny,indhmx,alength(0,0,2),value)
	call mx2indx(0,1,nxdim,nydim+1,nx,ny,indhmx,alength(0,0,1),value)

	end

c*******************************************************************

	subroutine rf49(nx,ny,nmax)

c reads header of unit 49

	implicit none

	integer nx,ny
	integer nmax

	include 'regflux.h'

	integer k

	read(49,*) nx,ny,nmax
	indhmx = nmax

	read(49,*) (iboxes(1,k),iboxes(2,k),k=1,indhmx)

	call indx2mx(1,1,nxdim,nydim,nx,ny,indhmx,hmed,value)
	call indx2mx(1,1,nxdim,nydim,nx,ny,indhmx,areav,value)
	call indx2mx(1,1,nxdim,nydim,nx,ny,indhmx,dxside,value)
	call indx2mx(1,1,nxdim,nydim,nx,ny,indhmx,dyside,value)
	call indx2mx(0,0,nxdim,nydim,nx,ny,indhmx,alength(0,0,2),value)
	call indx2mx(0,1,nxdim,nydim+1,nx,ny,indhmx,alength(0,0,1),value)

	end

c*******************************************************************

	subroutine wr49(it,nx,ny)

c writes data to unit 49

	implicit none

	integer it
	integer nx,ny

	include 'regflux.h'

	write(49,*) it

	call mx2indx(1,1,nxdim,nydim,nx,ny,indhmx,zlev,value)
	call mx2indx(0,0,nxdim,nydim,nx,ny,indhmx,fflux(0,0,2),value)
	call mx2indx(0,1,nxdim,nydim+1,nx,ny,indhmx,fflux(0,0,1),value)

	end

c*******************************************************************

	subroutine rd49(it,nx,ny,bstop)

c reads data from unit 49

	implicit none

	integer it
	integer nx,ny
	logical bstop

	include 'regflux.h'

	bstop = .false.

	read(49,*,end=99) it

	call indx2mx(1,1,nxdim,nydim,nx,ny,indhmx,zlev,value)
	call indx2mx(0,0,nxdim,nydim,nx,ny,indhmx,fflux(0,0,2),value)
	call indx2mx(0,1,nxdim,nydim+1,nx,ny,indhmx,fflux(0,0,1),value)

	return
   99	continue
	bstop = .true.
	return
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine mx2indx(nx0,ny0,nxd,nyd,nx,ny,nmax,source,target)

c copies 2D to 1D array

	implicit none

	include 'regflux.h'

	integer nx0,ny0
	integer nxd,nyd
	integer nx,ny
	integer nmax
	real source(nx0:nxd,ny0:nyd)
	real target(1)

	integer i,j,k

        do j=1,ny
          do i=1,nx
            target( gindsh(i,j) ) = source(i,j)
          end do
        end do

        write(49,*) (target(k),k=1,nmax)

	end

c*******************************************************************

	subroutine indx2mx(nx0,ny0,nxd,nyd,nx,ny,nmax,target,source)

c copies 1D to 2D array

	implicit none

	include 'regflux.h'

	integer nx0,ny0
	integer nxd,nyd
	integer nx,ny
	integer nmax
	real target(nx0:nxd,ny0:nyd)
	real source(1)

	integer i,j,k

        read(49,*) (source(k),k=1,nmax)

        do j=1,ny
          do i=1,nx
            target(i,j) = source( gindsh(i,j) )
          end do
        end do

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

        subroutine rdsind

c reads index from susanne

        implicit none

	include 'regflux.h'

        integer IG, KG, JE, INDHLD, INDHSE
        integer I,K
	integer index

	write(6,*) 'reading index field...'

        OPEN (60,FILE='GIND2D',status='old')
        read (60,'(20I6)') IG, KG, JE, INDHMX, INDHLD, INDHSE
        DO I=1,IG
          read (60,'(20I6)') (GINDH(I,K), K=1,KG)
        END DO
        CLOSE (60)

        do i=0,nydim
          do k=0,nxdim
            gindsh(k,i) = INDHLD
          end do
        end do

        do i=1,ig
          do k=1,kg
            gindsh(k,i) = GINDH(I,K)
          end do
        end do

	write(6,*) 'indhmx,maxind: ',indhmx,maxind

        if( indhmx .gt. maxind ) then
          write(6,*) 'indhmx,maxind: ',indhmx,maxind
          stop 'error stop rdsind: maxind'
        end if

c inverse index

	do i=1,indhmx
	  iboxes(1,i) = 0
	  iboxes(2,i) = 0
	end do

        do i=1,ig
          do k=1,kg
            index = gindsh(k,i)
	    iboxes(1,index) = k
	    iboxes(2,index) = i
          end do
        end do

        end

c*******************************************************************

