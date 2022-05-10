
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
!    Copyright (C) 2015  Erik Pascolo
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

c routines to create subsets for OMP
c
c contents :
c
c revision log :
c
c 01.09.2015	erp	code written from scratch
c 18.09.2015	ggu	code integrated
c 29.09.2015	ggu	new routine which is much faster
c 30.09.2015	ggu	new routine domain_clusterization_dummy()
c 05.10.2015	ggu	equilibrate subset filling
c 12.10.2015	ggu	changed VERS_7_3_3
c 16.11.2015	ggu	changed VERS_7_3_14
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c
c***********************************************************
c***********************************************************
c***********************************************************

!==================================================================
        module mod_subset
!==================================================================

        implicit none

        integer, save :: subset_num		!total number of subsets
        integer, save :: max_sharing_node  

        integer,save,allocatable,dimension(:)   :: subset_el !elems in subset
        integer,save,allocatable,dimension(:,:) :: indipendent_subset !elems
        integer,save,allocatable,dimension(:,:) :: nodes_map
        integer,save,allocatable,dimension(:,:) :: nodes_counter

!==================================================================
        contains
!==================================================================

      subroutine domain_clusterization

      use basin

      !call domain_clusterization_eric
      call domain_clusterization_ggu
      !call domain_clusterization_dummy

      end subroutine domain_clusterization

c***********************************************************

      subroutine domain_clusterization_eric
      
!$	use omp_lib
	use basin
    	
      implicit none
      
      integer :: ie,ik,i,j,k,n,m,a,b,c
      integer :: row_cluster
      integer :: numsubset,start,subset_lenght
      integer :: sum_subset_lenght
      integer :: max_subset_lenght
      integer :: max_l,min_l
      logical :: not_contained,stop_criterion
      integer,allocatable,dimension(:) :: subset
      integer, dimension(2) :: temp
      double precision :: timer
      
      !----------------------------------------------
      !! CREATE NODE COUNTER  
      !----------------------------------------------
      
      print *,"  ----------- STARTING DOMAIN CLUSTERING -----------"

      allocate(nodes_counter(2,nkn))
            
      nodes_counter = 0
     
      do ie=1,nel
         nodes_counter(1,nen3v(1,ie)) = nen3v(1,ie)
	 nodes_counter(2,nen3v(1,ie)) = nodes_counter(2,nen3v(1,ie))+1
	 nodes_counter(1,nen3v(2,ie)) = nen3v(2,ie)
	 nodes_counter(2,nen3v(2,ie)) = nodes_counter(2,nen3v(2,ie))+1
	 nodes_counter(1,nen3v(3,ie)) = nen3v(3,ie)
	 nodes_counter(2,nen3v(3,ie)) = nodes_counter(2,nen3v(3,ie))+1
      enddo
      
      max_sharing_node = MAXVAL(nodes_counter(2,:))
      
      !! sort nodes_counter in descending order

      DO j = 1,nkn-1
	DO k = j+1,nkn
	  IF(nodes_counter(2,j) .lt. nodes_counter(2,k)) THEN
	  temp = nodes_counter(:,k)
	  nodes_counter(:,k) =nodes_counter(:,j)
	  nodes_counter(:,j) = temp
	  ENDIF
	END DO 
      END DO 

      !----------------------------------------------
      !! END CREATE NODE COUNTER
      !----------------------------------------------
      
      !----------------------------------------------
      !! CREATE ELEMENTS MAPS
      !----------------------------------------------
      
      allocate(nodes_map(0:max_sharing_node,nkn))
      nodes_map = 0
      
      do i=1,nkn
	n = 1
	nodes_map(0,i) = i
	do j=1,nel
	  
	  if(nen3v(1,j) .eq. i .OR. nen3v(2,j) .eq. i .OR.
     +  	  nen3v(3,j) .eq. i) then
	    nodes_map(n,i) = j
	    n = n+1
	  endif
	  
	end do
      end do
      
      !----------------------------------------------
      !! END CREATE ELEMENTS MAPS
      !----------------------------------------------
      
      print *,"  ----------- DOMAIN CLUSTERING INFORMATION -----------"
      print *,"  NUM NODES = ",nkn," NUM ELEMENTS = ",nel
      print *,"  MAX LINK = ",max_sharing_node
      
      !----------------------------------------------
      !! GREEDY LOOP
      !----------------------------------------------
      
      allocate(subset(nel))
      
!$    timer = omp_get_wtime()
      
      subset = 0
      numsubset = 1
      max_subset_lenght = 0
      sum_subset_lenght = 0
      start = 1
      stop_criterion = .true.
      
      do while(stop_criterion)  
            
      call greedy_subset(start,numsubset,subset_lenght,subset,
     +                       nel,max_l,min_l)
      
      print *,"  SUBSET = ",numsubset," LENGTH = ",subset_lenght
      
      max_subset_lenght = MAX(max_subset_lenght,subset_lenght)
      sum_subset_lenght = sum_subset_lenght + subset_lenght
      
      stop_criterion = .false.
      
     	do a=1,nel
	  IF(subset(a) == 0) then
	      start = a
	      stop_criterion = .true.
	      numsubset = numsubset + 1
	      exit
	  endif
        end do 
            
      end do
      
      !----------------------------------------------
      !! END GREEDY LOOP
      !----------------------------------------------
     
!$     timer = omp_get_wtime() - timer
      
      !----------------------------------------------
      !! CREATE INDIPENDENT SUBSET MATRIX
      !----------------------------------------------

      subset_num = numsubset
      allocate(indipendent_subset(max_subset_lenght,subset_num))
      allocate(subset_el(subset_num))
      indipendent_subset = 0
  
     !! fill subset matrix loop     
      do i=1,subset_num
	  n = 1
	  do j = 1,nel
	      if(subset(j) == i) then
		indipendent_subset(n,i) = j
		n = n+1
	      endif
	    enddo
	    subset_el(i) = n-1
      enddo
      
      print *,"  NUM ITERATION GREEDY ALGO = ",numsubset
      print *,"  SUM SUBSET LENGTH = ",sum_subset_lenght
      print *,"  MAX LENGTH SUBSET = ",max_subset_lenght
      print *,"  TIME GREEDY = ",timer
      
      !----------------------------------------------
      ! controll if subset is indipendent
      !----------------------------------------------

      call check_subset(subset,max_subset_lenght,nel)
      
      print *,"  -----------------------------------------------------"
      
      deallocate(subset)
     
      end subroutine domain_clusterization_eric

c***********************************************************
      
      subroutine greedy_subset(start,numsubset,subset_lenght,subset
     +                           ,nel,max_l,min_l)
      
      implicit none
      
      integer,intent(in) :: nel,start,numsubset
      integer,intent(out) :: subset_lenght,max_l,min_l
      integer :: i,j,k,jend,jel,check_indipendency_global
      integer :: check_indipendency_local
      integer,dimension(:),intent(inout) :: subset

      max_l = 0
      min_l = 1000
      subset(start) = numsubset
      subset_lenght = 1
      check_indipendency_global = 0
      
      do i = 1,nel 
	
	! check if the element is just added
	if( subset(i) == 0 ) then
	
	 check_indipendency_global = 0
	  
	   ! loop on just assigned element
!$OMP PARALLEL DO PRIVATE(j)  
!$OMP& FIRSTPRIVATE(numsubset,i) DEFAULT(NONE)
!$OMP& SHARED(subset,nel) SCHEDULE(GUIDED)
!$OMP& REDUCTION(+:check_indipendency_global)
	   do j=1,nel
	   
            ! check all elements that have the same numsubset 
	      if(subset(j) == numsubset)then
	    
		! check indipendency with other subset element
		if(indipendent_element(i,j) .eqv. .false.) then
		  check_indipendency_global=check_indipendency_global+1
		endif
	      
	      endif
	   
	  enddo
!$OMP END PARALLEL DO

	  ! if all elements in subset is indipendent add el i
	  if(check_indipendency_global == 0) then
	    subset(i) = numsubset
	    subset_lenght = subset_lenght + 1
	   
	  endif
	  
	endif
            
      enddo
      
      end subroutine greedy_subset

c***********************************************************

      logical function indipendent_element(i,j)
      
	use basin, only : nen3v

      implicit none
      
      integer,intent(in) :: i,j
      integer :: a,b,c
      !logical :: indipendent_element
      
      indipendent_element = .true.
      	   
 	      do b = 1,3
 	       do c = 1,3
		if(nen3v(b,i) .eq. nen3v(c,j)) then
 		indipendent_element = .false.
 		endif
 	      enddo
 	    enddo
      
      end function indipendent_element
      
c***********************************************************

      subroutine check_subset(elem,max_subset_el,nel)
      
      implicit none
      integer, intent(in) :: nel,max_subset_el
      integer,intent(inout),dimension(:) :: elem
      integer :: ie,ik,i,j,k,n,m
      
      do i=1,nel
	 IF(elem(i) .eq. 0) then
	 print *,"ERROR ",i," is not set"
	 endif
      enddo
      
      do i=1,subset_num

	do j=1,subset_el(i)
	  
	  ie = indipendent_subset(j,i)
	  
	  do k=1,subset_el(i)
	  
	  if( k .ne. j ) then
	    
	    ik = indipendent_subset(k,i)
 	    if (indipendent_element(ie,ik) .neqv. .true.) then
	       print *,"  ERROR SUBSET NOT INDEPENDENT",i,ie,ik
	       stop
	    endif
	    
	   endif
	  
	  enddo
	  
	enddo
      enddo

      print *,"  ALL SUBSETS ARE INDIPENDENT"
      
      end subroutine check_subset

c***********************************************************
c***********************************************************
c***********************************************************

	subroutine domain_clusterization_ggu

	use basin

	implicit none

        !integer, save :: subset_num		!total number of subsets
        !integer,save,allocatable,dimension(:)   :: subset_el !elems in subset
        !integer,save,allocatable,dimension(:,:) :: indipendent_subset !elems

	integer ie
	integer ncol,nmax,ncolor,naver
	integer ic,n,ntot
	integer colork(nkn)
	integer color(nel)

	nmax = 0
	color = 0
	colork = 0

        print *,"  ----------- STARTING DOMAIN CLUSTERING -----------"

!	-------------------------------------------
!	color elements
!	-------------------------------------------

	call compute_color(ncolor,nmax)
	naver = nint(nel/float(ncolor))
	write(6,*) 'ncolor=',ncolor,' average=',naver,' nmax=',nmax

!	-------------------------------------------
!	set parameters and allocate arrays
!	-------------------------------------------

	subset_num = ncolor
	allocate(subset_el(subset_num))
	subset_el = 0

	call compute_subset_filling(nmax)	!sets subset_el

!	-------------------------------------------
!	adjust colors
!	-------------------------------------------

	call adjust_color

!	-------------------------------------------
!	check colors
!	-------------------------------------------

	call check_color

!	-------------------------------------------
!	transfer complete information to arrays
!	-------------------------------------------

	call compute_subset_filling(nmax)	!sets subset_el
	allocate(indipendent_subset(nmax,subset_num))
	subset_el = 0
	indipendent_subset = 0

	do ie=1,nel
	  ic = color(ie)
	  n = subset_el(ic) + 1
	  if( n > nmax ) then
	    write(6,*) 'n,nmax: ',n,nmax
	    stop 'error stop domain_clusterization_ggu: intern (1)'
	  end if
	  indipendent_subset(n,ic) = ie
	  subset_el(ic) = n
	end do

!	-------------------------------------------
!	write to terminal
!	-------------------------------------------

        print *,"  ---------- DOMAIN CLUSTERING INFORMATION ----------"
        print *,"  NUM NODES = ",nkn," NUM ELEMENTS = ",nel
        print *,"  NUMBER OF SUBSETS = ",subset_num
        ntot = 0
        do ic=1,subset_num
  	  n = subset_el(ic)
	  ntot = ntot + n
          print *,"  SUBSET = ",ic," LENGTH = ",n
        end do
        print *,"  SUM SUBSET LENGTH = ",ntot

	if( ntot /= nel ) then
	  write(6,*) 'ntot,nel: ',ntot,nel
	  write(6,*) 'error in total length of subsets'
	  stop 'error stop domain_clusterization_ggu: internal error (1)'
	end if

!	-------------------------------------------
!	end of routine
!	-------------------------------------------

	contains

!***********************

	subroutine compute_ic_minmax(icmin,icmax)

	integer icmin,icmax

	integer ic,n
	integer ncmin,ncmax

	ncmin = nel
	ncmax = 0
	icmin = 0
	icmax = 0

	do ic=1,subset_num
	  n = subset_el(ic)
	  if( n > ncmax ) then
	    ncmax = n
	    icmax = ic
	  end if
	  if( n < ncmin ) then
	    ncmin = n
	    icmin = ic
	  end if
	end do

	end subroutine compute_ic_minmax

!***********************

	subroutine compute_subset_filling(nmax)

	integer nmax

	integer ie,ic,n

	nmax = 0
	subset_el = 0

	do ie=1,nel
	  ic = color(ie)
	  n = subset_el(ic) + 1
	  nmax = max(nmax,n)
	  subset_el(ic) = n
	end do

	end subroutine compute_subset_filling

!***********************

	subroutine adjust_color

	logical bdebug
	integer ic
	integer icmin,icmax,ncolor
	integer ncmin,ncmax
	integer naver,nemax
	integer ie_start

	bdebug = .true.
	bdebug = .false.
	ie_start = 0
	ncolor = subset_num
	naver = nint(nel/float(ncolor))
	icmax = 1
	icmin = ncolor

	write(6,'(20i6)') (subset_el(ic),ic=1,ncolor)

	do

	  call compute_ic_minmax(icmin,icmax)
	  ncmin = subset_el(icmin)
	  ncmax = subset_el(icmax)
	  nemax = min(ncmax-naver,naver-ncmin)
	  if( bdebug ) write(6,*) 'exchanging colors: ',icmin,icmax
	  if( bdebug ) write(6,*) 'aver and max changes: ',naver,nemax

	  call apply_one_color(icmax,icmin,ie_start,nemax,ncol)

	  call compute_subset_filling(nmax)	!sets subset_el

	  write(6,'(20i6)') (subset_el(ic),ic=1,ncolor)
	  if( bdebug ) write(6,*) 'changes done and max: ',ncol,nmax

	  if( ncol == 0 ) exit

	end do

	write(6,*) 'difference in filling: ',ncmin,ncmax,ncmax-ncmin
	end subroutine adjust_color

!***********************

	subroutine compute_color(ncolor,nmax)

	integer ncolor,nmax

	integer ncol,ic_new,ic_old
	integer ie_start,nemax

	nmax = 0
	ic_new = 0
	ic_old = 0
	ie_start = 0
	nemax = nel

	color = 0
	colork = 0

	do 
	  ic_new = ic_new + 1
	  if( ic_new > 31 ) stop 'error stop compute_color: no more bytes'
	  call apply_one_color(ic_old,ic_new,ie_start,nemax,ncol)
	  nmax = max(nmax,ncol)
	  if( ncol == 0 ) exit		!everything colored
	end do

	ncolor = ic_new - 1

	end subroutine compute_color

!***********************

	subroutine apply_one_color(ic_old,ic_new,ie_start,nemax,ncol)

	integer ic_old,ic_new	!colors to use
	integer ie_start	!starting point for ie
	integer nemax		!at most color this number of elements
	integer ncol		!how many elements colored

	integer n,ie,iestride

	ncol = 0
	n = 0
	ie = ie_start
	iestride = 1

	do
	  n = n + 1
	  if( n > nel ) exit
	  if( ncol >= nemax ) exit
	  ie = ie + iestride
	  if( ie > nel ) ie = mod(ie,nel)
	  if( may_color(ie,ic_old,ic_new) ) then
	    ncol = ncol + 1
	    call color_elem(ie,ic_new)
	  end if
	end do
	
	ie_start = ie

	end subroutine apply_one_color

!***********************

	function may_color(ie,ic_from,ic_to)

	logical may_color
	integer ic_from,ic_to

	integer ie,ires
	integer ii,k

	may_color = .false.

	if( color(ie) /= ic_from ) return
	do ii=1,3
	  k = nen3v(ii,ie)
	  ires = ibits(colork(k),ic_to,1)
	  !if( colork(k) /= 0 ) return
	  if( ires /= 0 ) return
	end do

	may_color = .true.

	end function may_color

!***********************

	subroutine color_elem(ie,ic_new)

	integer ie,ic_new

	integer ii,k,ic_old
	integer col

	ic_old = color(ie)
	color(ie) = ic_new
	do ii=1,3
	  k = nen3v(ii,ie)
	  col = colork(k)
	  col = ibclr(col,ic_old)
	  col = ibset(col,ic_new)
	  colork(k) = col
	end do

	end subroutine color_elem

!***********************

	subroutine check_color

	integer ie,ii,k
	integer ic,col,ires
	logical berror

	berror = .false.

	do ie=1,nel
	  ic = color(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    col = colork(k)
	    ires = ibits(col,ic,1)
	    if( ires == 0 ) then
	      write(6,*) 'no such color: ',ie,ii,k,ic
	      berror = .true.
	    end if
	    col = ibclr(col,ic)
	    colork(k) = col
	  end do
	end do

	do k=1,nkn
	  if( colork(k) /= 0 ) then
	    write(6,*) 'still color on node: ',k,colork(k)
	    berror = .true.
	  end if
	end do

	if( berror ) then
	  stop 'error stop check_color: internal error (1)'
	else
	  write(6,*) 'check_color passed...'
	end if

	end subroutine check_color

!***********************

	end subroutine domain_clusterization_ggu

c***********************************************************

	subroutine domain_clusterization_dummy

	use basin

	implicit none

	integer nmax,ie,ic
	integer ntot,n

!	-------------------------------------------
!	set parameters and allocate arrays
!	-------------------------------------------

	nmax = nel
	subset_num = 1
	allocate(subset_el(subset_num))
	allocate(indipendent_subset(nmax,subset_num))
	subset_el = 0
	indipendent_subset = 0

!	-------------------------------------------
!	transfer information to arrays
!	-------------------------------------------

	ic = 1
	do ie=1,nel
	  indipendent_subset(ie,ic) = ie
	end do
	subset_el(ic) = nel

!	-------------------------------------------
!	write to terminal
!	-------------------------------------------

        print *,"  ---------- DOMAIN CLUSTERING INFORMATION ----------"
        print *,"  NUM NODES = ",nkn," NUM ELEMENTS = ",nel
        print *,"  NUMBER OF SUBSETS = ",subset_num
        ntot = 0
        do ic=1,subset_num
  	  n = subset_el(ic)
	  ntot = ntot + n
          print *,"  SUBSET = ",ic," LENGTH = ",n
        end do
        print *,"  SUM SUBSET LENGTH = ",ntot

!	-------------------------------------------
!	end of routine
!	-------------------------------------------

	end subroutine domain_clusterization_dummy

!==================================================================
        end module mod_subset
!==================================================================


