!-----------------------------------------------------------------
module util_part
!-----------------------------------------------------------------

  use basin
  use geom
  use levels
  use links
  use initialize_par

  contains

!-----------------------------------------------------------------

  subroutine domain_setup

    implicit none

    !local variables
    integer inlk

    inlk = 3*neldi+2*nkndi

    call mod_geom_init(nkndi,neldi,ngr)

    call mklenk(inlk,nkndi,neldi,nen3v,ilinkv,lenkv)

    call mklink(nkndi,ilinkv,lenkv,linkv)

    call set_ilhv

    call set_ilhkv

    return

  end subroutine domain_setup

!------------------------------------------------------------------

  subroutine refinement_part(kepart,knparts)

    implicit none

    !dummy arguments
    integer                             :: kepart(neldi),knparts

    !local variables
    integer, allocatable                :: ielems_order(:)
    integer,allocatable                 :: ivertices(:)
    integer, allocatable                :: irelocate(:)
    integer, allocatable                :: imynodes(:,:)
    integer                             :: irank,ielem,ienum,ihave
    integer                             :: ietmp,iidx0,iidx1,iidx
    integer                             :: ik0,ik1,ie0,ie1
    integer                             :: icounter,itomove
    integer                             :: ielems(12)
    logical                             :: ghave_elems(12)
    logical                             :: gbound
    
    !loop control
    integer                             :: jnidx,jeidx,jrank,jke

    allocate(ielems_order(neldi))
    allocate(ivertices(0:knparts))
   
    ivertices=0
    ielem=0
    do jrank=0,knparts-1
      ienum=0
      do jeidx=1,neldi
        if((kepart(jeidx)-1).eq.jrank) then
          ielem = ielem + 1
          ienum = ienum +1
          ielems_order(ielem) = jeidx
        end if
      end do
      ivertices(jrank+1) = ivertices(jrank) + ienum
    end do

    write(6,*) (ivertices(jrank+1)-ivertices(jrank),jrank=0,knparts-1)

    allocate(irelocate(4))
    allocate(imynodes(knparts,nkndi))

    do jrank=1,knparts
    itomove=-1
    do while(itomove.ne.0)
      irelocate=-1
      itomove=0
master: do jeidx=ivertices(jrank-1)+1,ivertices(jrank)
	  do jke=1,3
	    ik0=nen3v(jke,ielems_order(jeidx))
            iidx0=ilinkv(ik0)+1
            iidx1=ilinkv(ik0+1)
	    gbound=.false.
            if( lenkv(iidx1) .le. 0 ) then
	      gbound=.true.
	      iidx1 = iidx1 - 1
	    end if

	    ienum=0
            ielems=0
	    ghave_elems=.false.
            do iidx=iidx0,iidx1
	      ienum=ienum+1
              ielems(ienum) = lenkv(iidx)
	      do ihave=ivertices(jrank-1)+1,ivertices(jrank)
	        if(ielems_order(ihave).eq.ielems(ienum)) then
		  ghave_elems(ienum)=.true.
		end if
	      end do
	    end do

	    icounter=0
	    do iidx=1,ienum
	      if((ghave_elems(iidx))) then
		icounter=icounter+1
              end if
	    end do

	    if(icounter.eq.(ienum-1)) then

              itomove=itomove+1
              irelocate(3)=jrank-1
	      
     	      do iidx=1,ienum
                if(.not.ghave_elems(iidx)) then
		  irelocate(2)=ielems(iidx)
           	  exit
     	        end if
   	      end do

slave:	      do iidx0=1,knparts 
        	do iidx1=ivertices(iidx0-1)+1,ivertices(iidx0)
          	  if(ielems_order(iidx1).eq.irelocate(2)) then
            	    irelocate(4)=iidx0-1
            	    irelocate(1)=iidx1
            	    exit slave
          	  end if
        	end do
      	      end do slave

       	      if(irelocate(3).eq.(irelocate(4)+1)) then
	  	ietmp=ielems_order(ivertices(irelocate(4)+1))
	  	ielems_order(ivertices(irelocate(4)+1))=ielems_order(irelocate(1))
	  	ielems_order(irelocate(1))=ietmp
	  	if((irelocate(4)+1).gt.0) then
	    	  ivertices(irelocate(4)+1)=ivertices(irelocate(4)+1)-1
	  	end if
	  	cycle master
	      else if(irelocate(3).eq.(irelocate(4)-1)) then
	  	if(ivertices(irelocate(3)+1)+1.eq.irelocate(1)) then
	    	  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
		  cycle master
		end if
		ietmp=ielems_order(ivertices(irelocate(3)+1)+1)
		ielems_order(ivertices(irelocate(3)+1)+1)=ielems_order(irelocate(1))
		ielems_order(irelocate(1))=ietmp
		ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
		cycle master
	      else if(irelocate(3).lt.irelocate(4)) then
		ietmp=ielems_order(irelocate(1))
		ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4))+1)
	      	do iidx0=irelocate(4),irelocate(3)+2,-1
		  ielems_order(ivertices(iidx0)+1)=ielems_order(ivertices(iidx0-1)+1)
		end do
	      	do iidx0=irelocate(4),irelocate(3)+2,-1
		  ivertices(iidx0)=ivertices(iidx0)+1
		end do
		ielems_order(ivertices(irelocate(3)+1)+1)=ietmp
		ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
		cycle master
	      else if(irelocate(3).gt.irelocate(4)) then
		ietmp=ielems_order(irelocate(1))
		ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4)+1))
	      	do iidx0=irelocate(4)+1,irelocate(3)-1
		  ielems_order(ivertices(iidx0))=ielems_order(ivertices(iidx0+1))
		end do
	      	do iidx0=irelocate(4)+1,irelocate(3)-1
		  ivertices(iidx0)=ivertices(iidx0)-1
		end do
		ielems_order(ivertices(irelocate(3)))=ietmp
		ivertices(irelocate(3))=ivertices(irelocate(3))-1
		cycle master
	      end if

	    end if
	    
	  end do
	end do master
      end do 
    end do !close external loop

    itomove=-1

      irelocate=-1
      itomove=0
      do jrank=1, knparts
        ielem=1
	do jeidx=ivertices(jrank-1)+1,ivertices(jrank)
	  do jke=1,3
	    ik0=nen3v(jke,ielems_order(jeidx))
	    if(lenkv(ilinkv(ik0+1)).ne.0) cycle
	    do iidx0=1,ielem
	        if(imynodes(jrank,iidx0) .eq. ik0) then
		  exit
	        end if
	    end do
 	    if (iidx0 .gt. ielem) then
	      imynodes(jrank,ielem)=ik0
	      ielem=ielem+1
	    end if
	  end do
        end do

	do jnidx=1,ielem
	  ik0=imynodes(jrank,jnidx)
	  ienum=ilinkv(ik0+1)-ilinkv(ik0)
master4:  do jke=1,ienum
	    ik1=linkv(ilinkv(ik0)+jke)
	    do iidx0=1,ielem
	      if(ik1.eq.imynodes(jrank,iidx0)) exit
	    end do
	    if(iidx0 .gt. ielem) cycle
	    if(jke.eq.1) then
	      ie0= lenkv(ilinkv(ik0)+jke)
	      ie1= lenkv(ilinkv(ik0)+ienum)
	      do iidx0=ivertices(jrank-1)+1,ivertices(jrank)
		if(ie0.eq.ielems_order(iidx0)) exit
		if(ie1.eq.ielems_order(iidx0)) exit
	      end do
	      if (iidx0.gt.ivertices(jrank)) then
		irelocate(3)=jrank-1
		if(ie0.ne.0) then
		  irelocate(2)=ie0
		else if(ie1.ne.0) then
		  irelocate(2)=ie1
		else
		  write(6,*)'error not_elems',ie0,ie1
		end if

slave4:         do iidx0=1,knparts
                  do iidx1=ivertices(iidx0-1)+1,ivertices(iidx0)
                    if(ielems_order(iidx1).eq.irelocate(2)) then
                      irelocate(4)=iidx0-1
                      irelocate(1)=iidx1
                      exit slave4
                    end if
                  end do
                end do slave4

                if(irelocate(3).eq.(irelocate(4)+1)) then
                  ietmp=ielems_order(ivertices(irelocate(4)+1))
                  ielems_order(ivertices(irelocate(4)+1))=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  if((irelocate(4)+1).gt.0) then
                  ivertices(irelocate(4)+1)=ivertices(irelocate(4)+1)-1
                end if
                exit master4
                else if(irelocate(3).eq.(irelocate(4)-1)) then
                  if(ivertices(irelocate(3)+1)+1.eq.irelocate(1)) then
                    ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                    exit master4
                  end if
                  ietmp=ielems_order(ivertices(irelocate(3)+1)+1)
                  ielems_order(ivertices(irelocate(3)+1)+1)=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).lt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4))+1)
                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ielems_order(ivertices(iidx0)+1)=ielems_order(ivertices(iidx0-1)+1)
                  end do

                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ivertices(iidx0)=ivertices(iidx0)+1
                  end do

                  ielems_order(ivertices(irelocate(3)+1)+1)=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).gt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4)+1))
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ielems_order(ivertices(iidx0))=ielems_order(ivertices(iidx0+1))
                  end do
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ivertices(iidx0)=ivertices(iidx0)-1
                  end do
                  ielems_order(ivertices(irelocate(3)))=ietmp
                  ivertices(irelocate(3))=ivertices(irelocate(3))-1
                  exit master4
                end if



	      end if
	    else if(jke.eq.ienum) then
	      ie0= lenkv(ilinkv(ik0)+ienum-1)
	      ie1= lenkv(ilinkv(ik0)+ienum)
	      do iidx0=ivertices(jrank-1)+1,ivertices(jrank)
		if(ie0.eq.ielems_order(iidx0)) exit
		if(ie1.eq.ielems_order(iidx0)) exit
	      end do
	      if (iidx0.gt.ivertices(jrank)) then
		irelocate(3)=jrank-1
		if(ie0.ne.0) then
		  irelocate(2)=ie0
		else if(ie1.ne.0) then
		  irelocate(2)=ie1
		else
		  write(6,*)'error not_elems',ie0,ie1
		end if

                do iidx0=1,knparts
                  do iidx1=ivertices(iidx0-1)+1,ivertices(iidx0)
                    if(ielems_order(iidx1).eq.irelocate(2)) then
                      irelocate(4)=iidx0-1
                      irelocate(1)=iidx1
                      exit
                    end if
                  end do
                end do

                if(irelocate(3).eq.(irelocate(4)+1)) then
                  ietmp=ielems_order(ivertices(irelocate(4)+1))
                  ielems_order(ivertices(irelocate(4)+1))=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  if((irelocate(4)+1).gt.0) then
                  ivertices(irelocate(4)+1)=ivertices(irelocate(4)+1)-1
                end if
                exit master4
                else if(irelocate(3).eq.(irelocate(4)-1)) then
                  if(ivertices(irelocate(3)+1)+1.eq.irelocate(1)) then
                    ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                    exit master4
                  end if
                  ietmp=ielems_order(ivertices(irelocate(3)+1)+1)
                  ielems_order(ivertices(irelocate(3)+1)+1)=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).lt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4))+1)
                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ielems_order(ivertices(iidx0)+1)=ielems_order(ivertices(iidx0-1)+1)
                  end do

                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ivertices(iidx0)=ivertices(iidx0)+1
                  end do

                  ielems_order(ivertices(irelocate(3)+1)+1)=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).gt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4)+1))
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ielems_order(ivertices(iidx0))=ielems_order(ivertices(iidx0+1))
                  end do
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ivertices(iidx0)=ivertices(iidx0)-1
                  end do
                  ielems_order(ivertices(irelocate(3)))=ietmp
                  ivertices(irelocate(3))=ivertices(irelocate(3))-1
                  exit master4
                end if


	      end if
	    else
	      ie0= lenkv(ilinkv(ik0)+jke-1)
	      ie1= lenkv(ilinkv(ik0)+jke)
	      do iidx0=ivertices(jrank-1)+1,ivertices(jrank)
		if(ie0.eq.ielems_order(iidx0)) exit
		if(ie1.eq.ielems_order(iidx0)) exit
	      end do
	      if (iidx0.gt.ivertices(jrank)) then

		irelocate(3)=jrank-1
		if(ie0.ne.0) then
		  irelocate(2)=ie0
		else if(ie1.ne.0) then
		  irelocate(2)=ie1
		else
		  write(6,*)'error not_elems',ie0,ie1
		end if

                do iidx0=1,knparts
                  do iidx1=ivertices(iidx0-1)+1,ivertices(iidx0)
                    if(ielems_order(iidx1).eq.irelocate(2)) then
                      irelocate(4)=iidx0-1
                      irelocate(1)=iidx1
                      exit
                    end if
                  end do
                end do

                if(irelocate(3).eq.(irelocate(4)+1)) then
                  ietmp=ielems_order(ivertices(irelocate(4)+1))
                  ielems_order(ivertices(irelocate(4)+1))=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  if((irelocate(4)+1).gt.0) then
                  ivertices(irelocate(4)+1)=ivertices(irelocate(4)+1)-1
                end if
                exit master4
                else if(irelocate(3).eq.(irelocate(4)-1)) then
                  if(ivertices(irelocate(3)+1)+1.eq.irelocate(1)) then
                    ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                    exit master4
                  end if
                  ietmp=ielems_order(ivertices(irelocate(3)+1)+1)
                  ielems_order(ivertices(irelocate(3)+1)+1)=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).lt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4))+1)
                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ielems_order(ivertices(iidx0)+1)=ielems_order(ivertices(iidx0-1)+1)
                  end do

                  do iidx0=irelocate(4),irelocate(3)+2,-1
                    ivertices(iidx0)=ivertices(iidx0)+1
                  end do

                  ielems_order(ivertices(irelocate(3)+1)+1)=ietmp
                  ivertices(irelocate(3)+1)=ivertices(irelocate(3)+1)+1
                  exit master4
                else if(irelocate(3).gt.irelocate(4)) then
                  ietmp=ielems_order(irelocate(1))
                  ielems_order(irelocate(1))=ielems_order(ivertices(irelocate(4)+1))
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ielems_order(ivertices(iidx0))=ielems_order(ivertices(iidx0+1))
                  end do
                  do iidx0=irelocate(4)+1,irelocate(3)-1
                    ivertices(iidx0)=ivertices(iidx0)-1
                  end do
                  ielems_order(ivertices(irelocate(3)))=ietmp
                  ivertices(irelocate(3))=ivertices(irelocate(3))-1
                  exit master4
                end if



	      end if
	    end if
	  end do master4
      end do

    end do !close external loop

    itomove=-1
    do while(itomove.ne.0)

      irelocate=-1
      itomove=0

      call mv_elems(ielems_order,ivertices,irelocate,knparts,itomove)

    end do

   
   write(6,*) '----------------- End Refinement Partitioning -------------------'
        
    write(6,*) (ivertices(jrank+1)-ivertices(jrank),jrank=0,knparts-1)

    irank=-1
    do jeidx=1,neldi
      do jrank=1,knparts
	if(jeidx.le.ivertices(jrank)) then
	  irank=jrank-1
	  exit
	end if 
      end do
      kepart(ielems_order(jeidx)) = irank
    end do

    return

  end subroutine

!-----------------------------------------------------------------------------------

  subroutine mv_elems(kelems_order,kvertices,krelocate,knparts,ktomove)

    implicit none

    !dummy arguments
    integer, intent(inout)              :: kelems_order(neldi)
    integer, intent(inout)              :: kvertices(0:knparts)
    integer, intent(inout)              :: krelocate(1:4)
    integer, intent(in)                 :: knparts
    integer, intent(inout)              :: ktomove

    !loop control
    integer                             :: jrank,jeidx,jke

    !local variables
    integer                             :: ik0,iidx,iidx0,iidx1,ieidx,ietmp
    integer                             :: ihave,icounter,icounter2
    integer                             :: inohave,ielect,iplace
    integer                             :: ielems(1:12)
    integer                             :: inoelems(1:10)
    logical                             :: gbound
    logical                             :: ghave_elems(1:12)

      do jrank=1,knparts
	do jeidx=kvertices(jrank-1)+1,kvertices(jrank)
master2:	  do jke=1,3
	    ik0=nen3v(jke,kelems_order(jeidx))
            iidx0=ilinkv(ik0)+1
            iidx1=ilinkv(ik0+1)
	    gbound=.false.
            if( lenkv(iidx1) .le. 0 ) then
	      gbound=.true.
	      iidx1 = iidx1 - 1
	    end if

	    ieidx=0
            ielems=0
	    ghave_elems=.false.
            do iidx=iidx0,iidx1
	      ieidx=ieidx+1
              ielems(ieidx) = lenkv(iidx)
	      do ihave=kvertices(jrank-1)+1,kvertices(jrank)
	        if(kelems_order(ihave).eq.ielems(ieidx)) then
		  ghave_elems(ieidx)=.true.
		end if
	      end do
	    end do

	    icounter=0
            inoelems=0
            icounter2=0
            inohave=0
	    do iidx=1,ieidx
	      if((ghave_elems(iidx))) then
		icounter=icounter+1
	        inoelems(iidx)=0
                icounter2=0
              else
                icounter2=icounter2+1
                inoelems(iidx)=icounter2
              end if
	    end do
	    if(.not.gbound) then
              if(.not.ghave_elems(ieidx)) inoelems(1)=inoelems(ieidx)+1
	      if(.not.ghave_elems(1).and.ghave_elems(ieidx)) inohave=inohave+1
	    else
	      if(.not.ghave_elems(1).and.ghave_elems(ieidx)) inohave=inohave+1
	      if(ghave_elems(1).and.(ghave_elems(ieidx))) inohave=inohave+1
            end if

	    if (icounter .eq. ieidx) cycle

	    if(icounter.gt.1) then

   	    do iidx=1,ieidx-1
              if((ghave_elems(iidx)).and.(.not.ghave_elems(iidx+1))) then
       		inohave=inohave+1
     	      end if
   	    end do

   	    if(inohave.ge.2) then
     	      !write(6,*)'caso ihave'
     	      ielect=0
     	      do iidx=1,ieidx-1
       		if(inoelems(iidx).gt.0) then
	 	  if((ielect.gt.0).and.(ielect.le.inoelems(iidx)))cycle
	   
         	  if(ghave_elems(iidx+1)) then
	   	    ielect=inoelems(iidx)
	   	    iplace=iidx
         	  end if
       		end if
     	      end do
     	      if(.not.ghave_elems(ieidx).and.ghave_elems(1).and.inoelems(ieidx).lt.ielect) then
       		ielect=inoelems(ieidx)
       		iplace=ieidx
     	      end if



	      if(ielect.eq.1) then

	        iidx=iplace
         	ktomove=ktomove+1
         	krelocate(2)=ielems(iidx)
         	krelocate(3)=jrank-1

slave2:      		do iidx0=1,knparts 
        	  do iidx1=kvertices(iidx0-1)+1,kvertices(iidx0)
          	    if(kelems_order(iidx1).eq.krelocate(2)) then
            	      krelocate(4)=iidx0-1
            	      krelocate(1)=iidx1
            	      exit slave2
          	    end if
        	  end do
      		end do slave2
        	if(krelocate(3).eq.(krelocate(4)+1)) then
	  	  ietmp=kelems_order(kvertices(krelocate(4)+1))
	  	  kelems_order(kvertices(krelocate(4)+1))=kelems_order(krelocate(1))
	  	  kelems_order(krelocate(1))=ietmp
	  	  if((krelocate(4)+1).gt.0) then
	    	  kvertices(krelocate(4)+1)=kvertices(krelocate(4)+1)-1
	  	end if
	  	cycle master2
		else if(krelocate(3).eq.(krelocate(4)-1)) then
	  	  if(kvertices(krelocate(3)+1)+1.eq.krelocate(1)) then
	    	    kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
	    	    cycle master2
	  	  end if
	  	  ietmp=kelems_order(kvertices(krelocate(3)+1)+1)
	  	  kelems_order(kvertices(krelocate(3)+1)+1)=kelems_order(krelocate(1))
	  	  kelems_order(krelocate(1))=ietmp
	  	  kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
	  	  cycle master2
		else if(krelocate(3).lt.krelocate(4)) then
	  	  ietmp=kelems_order(krelocate(1))
	  	  kelems_order(krelocate(1))=kelems_order(kvertices(krelocate(4))+1)
      	  	  do iidx0=krelocate(4),krelocate(3)+2,-1
	    	    kelems_order(kvertices(iidx0)+1)=kelems_order(kvertices(iidx0-1)+1)
	  	  end do

      	  	  do iidx0=krelocate(4),krelocate(3)+2,-1
	    	    kvertices(iidx0)=kvertices(iidx0)+1
	  	  end do

	  	  kelems_order(kvertices(krelocate(3)+1)+1)=ietmp
	  	  kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
	  	  cycle master2
		else if(krelocate(3).gt.krelocate(4)) then
	  	  ietmp=kelems_order(krelocate(1))
	  	  kelems_order(krelocate(1))=kelems_order(kvertices(krelocate(4)+1))
      	  	  do iidx0=krelocate(4)+1,krelocate(3)-1
	    	    kelems_order(kvertices(iidx0))=kelems_order(kvertices(iidx0+1))
	  	  end do
      	  	  do iidx0=krelocate(4)+1,krelocate(3)-1
	    	    kvertices(iidx0)=kvertices(iidx0)-1
	  	  end do
		  kelems_order(kvertices(krelocate(3)))=ietmp
		  kvertices(krelocate(3))=kvertices(krelocate(3))-1
		  cycle master2
		end if

	      else
		iidx=iplace
         	ktomove=ktomove+1
            	krelocate(2)=ielems(iplace)
         	krelocate(4)=jrank-1


		if(iidx.eq.(ieidx-1)) then
		  iidx=iidx+1
		else if(iidx.lt.(ieidx-1).and.(.not.ghave_elems(iidx+2))) then
		  iidx=iidx+1
		else
		  if((iplace-ielect).gt.0) then
		    iidx=iplace-ielect
		  else
		    iidx=ieidx+iplace-ielect
		  end if	    
		end if

        	do iidx1=kvertices(jrank-1)+1,kvertices(jrank)
          	  if(kelems_order(iidx1).eq.ielems(iidx)) then
            	    krelocate(1)=iidx1
            	    exit
          	  end if
        	end do

slave3:      	do iidx0=1,knparts 
        	  do iidx1=kvertices(iidx0-1)+1,kvertices(iidx0)
          	    if(kelems_order(iidx1).eq.ielems(iplace)) then
            	      krelocate(3)=iidx0-1
            	      exit slave3
          	    end if
        	  end do
      		end do slave3
        	if(krelocate(3).eq.(krelocate(4)+1)) then
	  	  ietmp=kelems_order(kvertices(krelocate(4)+1))
	  	  kelems_order(kvertices(krelocate(4)+1))=kelems_order(krelocate(1))
	  	  kelems_order(krelocate(1))=ietmp
	  	  if((krelocate(4)+1).gt.0) then
	    	    kvertices(krelocate(4)+1)=kvertices(krelocate(4)+1)-1
	  	  end if
	  	  cycle master2
		else if(krelocate(3).eq.(krelocate(4)-1)) then
	  	  if(kvertices(krelocate(3)+1)+1.eq.krelocate(1)) then
	    	    kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
		    cycle master2
		  end if
		  ietmp=kelems_order(kvertices(krelocate(3)+1)+1)
		  kelems_order(kvertices(krelocate(3)+1)+1)=kelems_order(krelocate(1))
		  kelems_order(krelocate(1))=ietmp
		  kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
		  cycle master2
		else if(krelocate(3).lt.krelocate(4)) then
		  ietmp=kelems_order(krelocate(1))
		  kelems_order(krelocate(1))=kelems_order(kvertices(krelocate(4))+1)
	      	  do iidx0=krelocate(4),krelocate(3)+2,-1
		    kelems_order(kvertices(iidx0)+1)=kelems_order(kvertices(iidx0-1)+1)
		  end do
	      	  do iidx0=krelocate(4),krelocate(3)+2,-1
		    kvertices(iidx0)=kvertices(iidx0)+1
		  end do
		  kelems_order(kvertices(krelocate(3)+1)+1)=ietmp
		  kvertices(krelocate(3)+1)=kvertices(krelocate(3)+1)+1
		  cycle master2
		else if(krelocate(3).gt.krelocate(4)) then
		  ietmp=kelems_order(krelocate(1))
		  kelems_order(krelocate(1))=kelems_order(kvertices(krelocate(4)+1))
	      	  do iidx0=krelocate(4)+1,krelocate(3)-1
		    kelems_order(kvertices(iidx0))=kelems_order(kvertices(iidx0+1))
		  end do
	      	  do iidx0=krelocate(4)+1,krelocate(3)-1
		    kvertices(iidx0)=kvertices(iidx0)-1
		  end do
		  kelems_order(kvertices(krelocate(3)))=ietmp
		  kvertices(krelocate(3))=kvertices(krelocate(3))-1
		  cycle master2
		end if
	      end if

   	    end if

	    end if
	    
	  end do master2
	end do
      end do


  end subroutine

!--------------------------------------------------------------------------!
!-------------------  start remove_dups subroutine  -----------------------!
!--------------------------------------------------------------------------!
!##########################################################################!
!* This subroutine sort the elements of a vector and delete the duplicate *!
!* the elements must be integer positive except the 0                     *!
!##########################################################################!

  subroutine remove_dups(array,length,res,dimout)

    implicit none
    
    integer, dimension(:) :: array         ! The input
    !integer :: res(size(array))  ! The output
    integer :: length            ! The number of unique elements in output
    integer :: i, j, k, h, temp
    integer, optional :: dimout
    integer,allocatable,dimension(:) :: res

    if (present(dimout) .and. (.not. allocated(res))) then
      allocate(res(dimout))
    else if (.not. allocated(res)) then
       write(6,*)'error in remove_dups'
    end if

    if( size(array) .le. 0) then
      length = 0
      return
    end if

    do i=1,size(res)
    !do i=1,size(array)
      res(i) = -1
    end do

    k = 1
    res(1) = array(1)
    outer: do i=2,size(array)
       if(array(i) .eq. 0) cycle outer
       do j=1,k
          if ((res(j) == array(i)) .or. (array(i) .eq. 0)) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       res(k) = array(i)

       h=k
       do while (res(h) .lt. res(h-1))
         temp = res(h-1)
         res(h-1) = res(h)
         res(h) = temp
         if(h .eq. 2) exit
         h = h-1
       end do

    end do outer

    length = k

    return

  end subroutine remove_dups

!-----------------------------------------------------------------
end module util_part
!-----------------------------------------------------------------
