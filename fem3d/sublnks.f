c
c $Id: sublnks.f,v 1.6 2009-05-21 09:24:00 georg Exp $
c
c topological set up routines
c
c contents :
c
c subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)
c			 computes parameters about total numbers
c subroutine newlnk
c			 administrates new topological routines
c function winkk(k)
c			 total angle at node k
c subroutine setnar(nar)
c			 sets number of areas
c subroutine setwnk(wink)
c			 sets total angle at boundary nodes
c subroutine arper
c			 computes statistics about area of submerged zones
c
c revision log :
c
c 01.08.2003	ggu	created from sublnk.f
c 13.08.2003	ggu	new name set_link_info for newlnk
c 17.02.2006	ggu	do not all anymore zvbnds()
c 28.04.2009    ggu     links re-structured
c
c*****************************************************************

        subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

c computes parameters about total numbers
c
c all numbers are refered to actual elements and nodes inside system
c
c nnkn   total number of nodes
c nnel   total number of elements
c nnbn   total number of boundary nodes
c nnli   total number of links
c nnis   total number of islands
c nnod   total number of not allowed junctions (should be 0)
c
c formulas :
c
c nel = 2*nkn - nbn - 2 + 2*nis - nod = nbn - 2 + 2*nin + 2*nis - nod
c nin = nkn - nbn
c nli = ( 3*nel + nbn + nod ) / 2
c
c the following is always true:
c
c nel < 2*nkn
c nli = nel + (1/2)*nel + nbn < nel + nkn + nbn <= nel + 2*nkn   so
c nli < 2*nkn + nel

	use mod_geom
	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

        implicit none

c arguments
        integer nnkn,nnel,nnbn,nnli,nnis,nnod
c local
        integer k,ie,n,i,ne
	integer elems(maxlnk)
        logical bin
c statement functions
        logical iskbnd,iskins,iseins
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskins(k) = inodv(k).ne.-2
        iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0

        nnod=0
        do k=1,nkn
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
        do k=1,nkn
          if( iskbnd(k) ) nnbn=nnbn+1
          if( iskins(k) ) nnkn=nnkn+1
        end do

        nnel=0
        do ie=1,nel
          if( iseins(ie) ) nnel=nnel+1
        end do

        nnis=(nnel+nnbn-2*nnkn+2+nnod)/2
        nnli=(3*nnel+nnbn+nnod)/2

        end

c****************************************************************

c        subroutine newlnk
        subroutine set_link_info

c administrates new topological routines
c
c ... iwegv has already been set

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	include 'femtime.h'
c local
        character*80 nam,dir,file
        real wink
        integer nnkn,nnel,nnbn,nnli,nnis,nnod,nnar,nnnis,nnnel
c function
        integer iround,ifileo
c save
        integer n88
        save n88
        data n88 /0/

        if(n88.eq.0) call getinfo(n88)          !get unit number

        call setnod

        call nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)

        call setnar(nnar)	!number of areas
        call setwnk(wink)	!total angle at boundary nodes

        nnnis = iround( (wink-nnbn*180.)/360. ) + nnar

        nnnel=2*nnkn-nnbn+2*(nnnis-nnar)-nnod

        write(n88,'(a,i10,10i6)') 'newlnk: ',it
     +          ,nnkn,nnel,nnbn,nnli,nnis,nnod,nnar
     +          ,nnnis,nnnel-nnel,nel-nnel

        end

c****************************************************************

        function winkk(k)

c total angle at node k
c
c k       node to links
c
	use mod_geom
	use mod_geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

c arguments
	real winkk
        integer k
c local
        integer i,n,ie
	integer elems(maxlnk)
        double precision w
c functions
        integer ithis

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

c****************************************************************
c
        subroutine setnar(nar)
c
c sets number of areas
c
c nar   number of areas (return)
c
	use mod_geom
	use mod_geom_dynamic
	use basin

        implicit none
c
c arguments
        integer nar
c local
        integer i,ie,ieo,ien,n1,n2
        logical btest
c functions
        integer inext
c data
        data btest /.false./
c
        do ie=1,nel
          if(iwegv(ie).gt.0) iwegv(ie)=iwegv(ie)+3
        end do
c
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
c              write(6,*) ieo,iwegv(ieo),i,ien,iwegv(ien)
c              read(5,'(i10)') n
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
c
        if(btest) then
          write(6,*) nar,n1,n2
          write(6,*)
        end if

        do ie=1,nel
          if(iwegv(ie).lt.3) then
            write(6,*) 'we forgot something : ',ie,iwegv(ie)
            stop 'error stop setnar'
          end if
          iwegv(ie)=iwegv(ie)-3
        end do
c
        return
        end
c
c****************************************************************
c
        subroutine setwnk(wink)
c
c sets total angle at boundary nodes
c
c wink  total angle (return)
c
	use mod_geom_dynamic
	use evgeom
	use basin

        implicit none
c
c arguments
        real wink
c local
        integer ie,ii,k
        double precision w
c functions
        logical iskbnd
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
c
c sum angles
c
        w=0.
        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              if( iskbnd(k) ) w=w+ev(10+ii,ie)
            end do
          end if
	end do
c
        wink=w
c
        return
        end
c
c****************************************************************
c
        subroutine arper
c
c computes statistics about area of submerged zones
c
c ... iwegv has already been set
c
	use mod_geom_dynamic
	use evgeom
	use basin, only : nkn,nel,ngr,mbw

        implicit none
c
	include 'femtime.h'
c local
        real arin,arout,artot,area
        integer ie
c functions
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

c        write(88,'(i8,f12.3,2f12.2,e12.4)') it,arin,arout,artot

        end

c****************************************************************

