c
c $Id: mod_subset.f,v 1.4 2009-04-03 16:38:23 georg Exp $
c
c routines to create subset
c
c contents :
c
c revision log :
c
c 01.09.2015	eps	code written from scratch
c 18.09.2015	ggu	code integrated
c 29.09.2015	ggu	new routine which is much faster
c 30.09.2015	ggu	new routine domain_clusterization_dummy()
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
      
      start_el_loop: 
     +	do a=1,nel
	  IF(subset(a) == 0) then
	      start = a
	      stop_criterion = .true.
	      numsubset = numsubset + 1
	      exit start_el_loop 
	  endif
        end do start_el_loop 
            
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

	subroutine domain_clusterization_ggu

	use basin

	implicit none

        !integer, save :: subset_num		!total number of subsets
        !integer,save,allocatable,dimension(:)   :: subset_el !elems in subset
        !integer,save,allocatable,dimension(:,:) :: indipendent_subset !elems

	integer ie
	integer icolor,ncol,nmax,ncolor
	integer ic,n,ntot
	integer colork(nkn)
	integer color(nel)

	icolor = 0
	color = 0
	nmax = 0

        print *,"  ----------- STARTING DOMAIN CLUSTERING -----------"

!	-------------------------------------------
!	color elements
!	-------------------------------------------

	call compute_color(ncolor,nmax,color)
	write(6,*) 'ncolor: ',ncolor

!	-------------------------------------------
!	set parameters and allocate arrays
!	-------------------------------------------

	subset_num = ncolor
	allocate(subset_el(subset_num))
	allocate(indipendent_subset(nmax,subset_num))
	subset_el = 0
	indipendent_subset = 0

!	-------------------------------------------
!	transfer information to arrays
!	-------------------------------------------

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

!	-------------------------------------------
!	end of routine
!	-------------------------------------------

	contains

!***********************

	subroutine compute_color(ncolor,nmax,color)

	integer ncolor,nmax
	integer color(nel)

	integer ncol,icolor
	integer ie
	integer iestart,iestride,n

	nmax = 0
	icolor = 0
	color = 0

	iestart = 0
	iestride = 1		!must be prime
	ie = iestart

	do 
	  ncol = 0
	  icolor = icolor + 1
	  colork = 0
	  !do ie=1,nel
	  n = 0
	  do
	    n = n + 1
	    if( n > nel ) exit
	    ie = ie + iestride
	    if( ie > nel ) ie = mod(ie,nel)
	    if( may_color(ie) ) then
	      ncol = ncol + 1
	      call color_elem(ie,icolor)
	    end if
	  end do
	  nmax = max(nmax,ncol)
	  !write(6,*) icolor,ncol,nmax
	  if( ncol == 0 ) exit		!everything colored
	end do

	ncolor = icolor - 1

	end subroutine compute_color

!***********************

	function may_color(ie)

	logical may_color
	integer ie
	integer ii,k

	may_color = .false.

	if( color(ie) /= 0 ) return
	do ii=1,3
	  k = nen3v(ii,ie)
	  if( colork(k) /= 0 ) return
	end do

	may_color = .true.

	end function may_color

!***********************

	subroutine color_elem(ie,icolor)

	integer ie,icolor
	integer ii,k

	color(ie) = icolor
	do ii=1,3
	  k = nen3v(ii,ie)
	  colork(k) = 1
	end do

	end subroutine color_elem

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


