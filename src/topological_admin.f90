!
! $Id: sublnks.f,v 1.6 2009-05-21 09:24:00 georg Exp $
!
! topological set up routines
!
! contents :
!
! subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)
!			 computes parameters about total numbers
! subroutine newlnk
!			 administrates new topological routines
! function winkk(k)
!			 total angle at node k
! subroutine setnar(nar)
!			 sets number of areas
! subroutine setwnk(wink)
!			 sets total angle at boundary nodes
! subroutine arper
!			 computes statistics about area of submerged zones
!
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 13.08.2003	ggu	new name set_link_info for newlnk
! 17.02.2006	ggu	do not all anymore zvbnds()
! 28.04.2009    ggu     links re-structured
!
!*****************************************************************
!------------------------------------------------------------------
        module topological_admin
!------------------------------------------------------------------
        contains
!------------------------------------------------------------------

        subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

! computes parameters about total numbers
!
! all numbers are refered to actual elements and nodes inside system
!
! nnkn   total number of nodes
! nnel   total number of elements
! nnbn   total number of boundary nodes
! nnli   total number of links
! nnis   total number of islands
! nnod   total number of not allowed junctions (should be 0)
!
! formulas :
!
! nel = 2*nkn - nbn - 2 + 2*nis - nod = nbn - 2 + 2*nin + 2*nis - nod
! nin = nkn - nbn
! nli = ( 3*nel + nbn + nod ) / 2
!
! the following is always true:
!
! nel < 2*nkn
! nli = nel + (1/2)*nel + nbn < nel + nkn + nbn <= nel + 2*nkn   so
! nli < 2*nkn + nel

	use geom
	use geom_dynamic
	use basin, only : nkn,nel,ngr,mbw
	use shympi
        use lnku

        implicit none

! arguments
        integer nnkn,nnel,nnbn,nnli,nnis,nnod
! local
        integer k,ie,n,i,ne,ntot,j
	integer elems(maxlnk)
        logical bin
! statement functions
        logical iskbnd,iskins,iseins
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskins(k) = inodv(k).ne.-2
        iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0

        if(shympi_partition_on_elements()) then
          ntot = nkn_inner      !SHYMPI_ELEM - should be total nodes to use
        else
          ntot = nkn
        end if

        nnod=0

        do k=1,ntot
	  call get_elems_around(k,maxlnk,ne,elems)

          n=0
	  ie = elems(ne)
          bin=iseins(ie)
          do i=1,ne
	    ie = elems(i)
            if( iseins(ie) ) then
              if( .not. bin ) n=n+1
            end if
            bin=iseins(ie)
          end do

          if(n.gt.1) nnod=nnod+n-1
        end do

        nnbn=0
        nnkn=0
        do k=1,ntot
          if( iskbnd(k) ) nnbn=nnbn+1
          if( iskins(k) ) nnkn=nnkn+1
        end do

        nnel=0
        do ie=1,nel
          if( iseins(ie) ) nnel=nnel+1
        end do

        nnel = shympi_sum(nnel)
        nnbn = shympi_sum(nnbn)
        nnod = shympi_sum(nnod)
        nnkn = shympi_sum(nnkn)

        nnis=(nnel+nnbn-2*nnkn+2+nnod)/2
        nnli=(3*nnel+nnbn+nnod)/2

        end

!****************************************************************

!        subroutine newlnk
        subroutine set_link_info

! administrates new topological routines
!
! ... iwegv has already been set

	use geom
	use basin, only : nkn,nel,ngr,mbw,neldi
	use shympi
        use utility
        use defnames

        implicit none

	include 'femtime.h'
! local
        character*80 nam,dir,file
        double precision wink
        integer nnkn,nnel,nnbn,nnli,nnis,nnod,nnar,nnnis,nnnel
! save
        integer n88
        save n88
        data n88 /0/

        if(n88.eq.0) call getinfo(n88)          !get unit number

        !call setnod

        call nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

        call setnar(nnar)	!number of areas
        call setwnk(wink)	!total angle at boundary nodes

        if(bmpi) then
          nnar = shympi_min(nnar) !- n_threads
        end if

        nnnis = iround( (wink-nnbn*180.0d0)/360.0d0 ) + nnar

        nnnel=2*nnkn-nnbn+2*(nnnis-nnar)-nnod

        if(shympi_is_master()) then
          write(n88,'(a,i8,4i10,6i6)') 'newlnk: ',it,nnkn,nnel,nnbn,nnli,nnis,nnod,nnar    &
     &          ,nnnis,nnnel-nnel,neldi-nnel
	end if

        end

!****************************************************************

        function winkk(k)

! total angle at node k
!
! k       node to links
!
	use geom
	use geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw
        use lnku

        implicit none

! arguments
	double precision winkk
        integer k
! local
        integer i,n,ie
	integer elems(maxlnk)
        double precision w

	call get_elems_around(k,maxlnk,n,elems)

        w=0.
        do i=1,n
          ie=elems(i)
          if( iwegv(ie).eq.0 ) then
            w = w + ev(10+ithis(k,ie),ie)
          end if
        end do

        winkk=w

        end

!****************************************************************
!
        subroutine setnar(nar)
!
! sets number of areas
!
! nar   number of areas (return)
!
	use geom
	use geom_dynamic
	use basin
	use shympi
        use lnku

        implicit none
!
! arguments
        integer nar
! local
        integer i,ie,ieo,ien,n1,n2
        logical btest
! data
        data btest /.false./
!
        do ie=1,nel
          if(iwegv(ie).gt.0) iwegv(ie)=iwegv(ie)+3
        end do
!
        n1=0
        n2=0
        nar=0
        do ie=1,nel
          if(iwegv(ie).eq.0) then
            i=1
            ieo=ie
            do while(iwegv(ieo).lt.3)
              iwegv(ieo)=iwegv(ieo)+1
              i=mod(i,3)+1
              ien=ieltv(i,ieo)
              n1=n1+1
!              write(6,*) ieo,iwegv(ieo),i,ien,iwegv(ien)
!              read(5,'(i10)') n
              if(ien.gt.0) then !possible to enter
               if(iwegv(ien).lt.3) then !possible to enter
                if(iwegv(ien).eq.0.or.iwegv(ieo).eq.3) then !may enter
                  i=inext( nen3v(mod(i,3)+1,ieo) , ien )
                  ieo=ien
                  n2=n2+1
                end if
               end if
              end if
            end do
            nar=nar+1
          end if
        end do
!
        if(btest) then
          write(6,*) nar,n1,n2
          write(6,*)
        end if

        !if(bmpi) then
        !nar = shympi_min(nar) !- n_threads
        !call shympi_comment('shympi_sum(nar)')
        !end if

        do ie=1,nel
          if(iwegv(ie).lt.3) then
            write(6,*) 'we forgot something : ',ie,iwegv(ie)
            stop 'error stop setnar'
          end if
          iwegv(ie)=iwegv(ie)-3
        end do
!
        return
        end
!
!****************************************************************
!
        subroutine setwnk(wink)
!
! sets total angle at boundary nodes
!
! wink  total angle (return)
!
	use geom_dynamic
	use evgeom
	use basin
	use shympi

        implicit none
!
! arguments
        double precision wink
! local
        integer ie,ii,k
        double precision w
! functions
        logical iskbnd
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
!
! sum angles
!
        w=0.
        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              if( iskbnd(k) ) w=w+ev(10+ii,ie)
            end do
          end if
	end do

        wink=w

        wink = shympi_sum(wink)
        !call shympi_comment('shympi_sum(wink)')

        return
        end

!****************************************************************
!
        subroutine arper
!
! computes statistics about area of submerged zones
!
! ... iwegv has already been set
!
	use geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

        implicit none
!
	include 'femtime.h'
! local
        double precision arin,arout,artot,area
        integer ie
! functions
        logical isein
        isein(ie) = iwegv(ie).eq.0

        arin=0.
        arout=0.
        artot=0.

        do ie=1,nel
          area=ev(10,ie)
          if(isein(ie)) then
            arin=arin+area
          else
            arout=arout+area
          end if
          artot=artot+area
        end do

        arin=100.*arin/artot
        arout=100.*arout/artot
        artot=12.*artot

!        write(88,'(i8,f12.3,2f12.2,e12.4)') it,arin,arout,artot

        end

!****************************************************************

!------------------------------------------------------------------
        end module topological_admin
!------------------------------------------------------------------
