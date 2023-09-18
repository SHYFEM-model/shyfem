!
! $Id: subsrt.f,v 1.3 1998/01/22 16:49:57 georg Exp $
!
! sort routines
!
! contents :
!
! subroutine sortir(ipoint,value,ndim)	indirect sort of double precision array
! subroutine sortii(ipoint,ivalue,ndim)	indirect sort of interger array
! subroutine sort(n,ra)			heapsort
! subroutine isort(n,ra,index)		heapsort (indirect)
! function locate(n,ixx,index,ix)	locates index for ix in array ixx
!
!***********************************************************
!----------------------------------------------------------------------
        module sort_util
!----------------------------------------------------------------------
        contains
!----------------------------------------------------------------------

	subroutine sortir(ipoint,value,ndim)

! indirect sort of double precision array
!
! ipoint	pointer to values (must be already defined)
! value		values to be sorted
! ndim		total number of values to be sorted

	implicit none

	integer ndim
	integer ipoint(ndim)
	double precision value(ndim)

	integer i,k,n,n1

	do i=1,ndim-1
	  do k=ndim,i+1,-1
	    n=ipoint(k)
	    n1=ipoint(k-1)
	    if(value(n).lt.value(n1)) then
		ipoint(k)=n1
		ipoint(k-1)=n
	    end if
	  end do
	end do

	end

!***********************************************************

	subroutine sortii(ipoint,ivalue,ndim)
!
! indirect sort of interger array
!
! ipoint	pointer to values (must be already defined)
! ivalue	values to be sorted
! ndim		total number of values to be sorted
!
	implicit none

	integer ndim
	integer ipoint(ndim),ivalue(ndim)

	integer i,k,n,n1

	do i=1,ndim-1
	  do k=ndim,i+1,-1
	    n=ipoint(k)
	    n1=ipoint(k-1)
	    if(ivalue(n).lt.ivalue(n1)) then
		ipoint(k)=n1
		ipoint(k-1)=n
	    end if
	  end do
	end do

	end

!***********************************************************

      subroutine sort(n,ra)

! heapsort (double precision)

      implicit none

      integer n
      double precision ra(n)
!      integer ra(n)

      integer i,ir,j,l
      double precision rra
!      integer rra

      if(n.lt.2) return

      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      continue
	if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end

!***********************************************************

      subroutine sort_int(n,ra)

! heapsort (integer)

      implicit none

      integer n
!      double precision ra(n)
      integer ra(n)

      integer i,ir,j,l
!      double precision rra
      integer rra

      if(n.lt.2) return

      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      continue
	if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end

!******************************************************************

      subroutine isort(n,ra,index)

! heapsort

      implicit none

      integer n
      integer ra(n),index(n)

      integer i,ir,j,l,irra
      integer rra

      do i=1,n
	index(i)=i
      end do

      if(n.lt.2) return

      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
	  irra=index(l)
          rra=ra(irra)
        else
	  irra=index(ir)
          rra=ra(irra)
          index(ir)=index(1)
          ir=ir-1
          if(ir.eq.1)then
            index(1)=irra
            return
          endif
        endif
        i=l
        j=l+l
20      continue
	if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(index(j)).lt.ra(index(j+1)))j=j+1
          endif
          if(rra.lt.ra(index(j)))then
            index(i)=index(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        index(i)=irra
      go to 10
      end

!******************************************************************

      function locate(n,ixx,index,ix)

! locates index for ix in array ixx ( ixx(index(i)) is sorted )

      implicit none

      integer locate
      integer n
      integer ixx(n)
      integer index(n)
      integer ix

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(ix.gt.ixx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      if(ju.le.0) then                    !empty list
        locate=0
      else if(ix.eq.ixx(index(ju))) then  !found
        locate=index(ju)
      else                                !not found
        locate=0
      end if

      return
      end

!**********************************************************

!----------------------------------------------------------------------
        end module sort_util
!----------------------------------------------------------------------
