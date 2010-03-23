c
c $Id: sublnk.f,v 1.22 2003/08/12 16:33:23 georg Exp $
c
c topological set up routines
c
c contents :
c
c subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)
c			 computes parameters about total numbers
c subroutine setnod
c			 sets array inodv
c
c function kthis(i,ie)	 gets node at position i in ie
c function knext(k,ie) 	 gets node after k in ie
c function kbhnd(k,ie) 	 gets node before k in ie
c function ithis(k,ie) 	 gets i of k in ie
c function inext(k,ie) 	 gets i after k in ie
c function ibhnd(k,ie)	 gets i before k in ie
c
c subroutine pntfla(k,ipf,ipl)
c			 gets first and last element in element link (absolute)
c subroutine pntfl(k,ipf,ipl)
c			 gets first and last element in element link of node k
c
c subroutine ntlink(ilinkv,lenkv)
c			 sets up vector with element links and a pointer to it
c
c subroutine setkan(kantv)
c			 sets up vector kantv
c subroutine setdxy(dxv,dyv)
c			 sets up vector dxv and dyv
c
c subroutine setelt	 sets up vector ieltv
c subroutine checkelt	 checks integrity of ieltv
c function ielt(i,ie)	 gets i'th neigbour of element ie
c
c subroutine setlin(linkv)
c			 sets up vector linkv
c subroutine getlin(k,link,nlink)
c			 gets links for node k
c
c subroutine tstlnk
c			 sets up vector with element links and a pointer to it
c
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
c revised ...07.92 by ggu   $$lump  - lumping of matrix
c revised 05.08.92 by ggu   $$ibtyp3 - implementation of ibtyp=3
c 27.03.1998	ggu	eliminated /bnd/, /irv/
c 27.04.1998	ggu	$$NKNEL - do not call nknel in tstlnk
c 08.05.1998	ggu	new routine pntfla -> absolute element index
c 20.05.1998	ggu	hard coded unit 88 substituted (use ifileo)
c 14.07.1998	ggu	$$ibtyp4 - boundary type 4 integrated
c 21.08.1998    ggu     file sublnk splitted into lnk/dry
c 27.12.2001    ggu     use info unit number in newlnk
c 31.07.2003    ggu     eliminated setwnv and winv(1)
c 31.07.2003    ggu     comodity functions is_internal_node etc.
c 10.08.2003    ggu     do not use anymore -> sublnka/s/u and links
c
c*****************************************************************
c
        subroutine nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)
c
c computes parameters about total numbers
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
c
        implicit none
c
c arguments
        integer nnkn,nnel,nnbn,nnli,nnis,nnod
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer inodv(1), iwegv(1)
        integer ilinkv(1), lenkv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /inodv/inodv, /iwegv/iwegv
        common /ilinkv/ilinkv, /lenkv/lenkv
c local
        integer k,ie,n,ip,ip0,ip1
        logical bin
c statement functions
        logical iskbnd,iskins,iseins
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iskins(k) = inodv(k).ne.-2
        iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0
c
        nnod=0
        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
c
          n=0
          bin=iseins(lenkv(ip1))
          do ip=ip0,ip1
            if( iseins(lenkv(ip)) ) then
              if( .not. bin ) n=n+1
            end if
            bin=iseins(lenkv(ip))
          end do
c
          if(n.gt.1) nnod=nnod+n-1
        end do
c
        nnbn=0
        nnkn=0
        do k=1,nkn
          if( iskbnd(k) ) nnbn=nnbn+1
          if( iskins(k) ) nnkn=nnkn+1
        end do
c
        nnel=0
        do ie=1,nel
          if( iseins(ie) ) nnel=nnel+1
        end do
c
        nnis=(nnel+nnbn-2*nnkn+2+nnod)/2
        nnli=(3*nnel+nnbn+nnod)/2
c
        return
        end
c
c****************************************************************
c
        subroutine setnod
c
c sets array inodv
c
c inodv 
c	 0: internal node  
c	>0: open boundary node
c       -1: boundary node  
c	-2: out of system
c
c if open boundary node, inodv(k) is number of boundary (ggu 15.11.2001)
c
        implicit none
c
c parameter
	real winmax
	parameter(winmax=359.8)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer iwegv(1),nen3v(3,1),inodv(1)
        real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /iwegv/iwegv, /nen3v/nen3v, /inodv/inodv
        common /ev/ev
c local
        integer ie,ii,k,n
        integer ibc,ibtyp
c functions
        integer ipext
	integer itybnd,nkbnds,kbnds
c equivalence
        real winkv(1)
        equivalence(inodv(1),winkv(1))
c
c initialize array to hold angles
c
        do k=1,nkn
          winkv(k)=0.
        end do
c
c sum angles
c
        do ie=1,nel
          if(iwegv(ie).eq.0) then !element is in system
            do ii=1,3
              k=nen3v(ii,ie)
              winkv(k)=winkv(k)+ev(10+ii,ie)
            end do
          end if
	end do
c
c set up inodv
c
        do k=1,nkn
          if(winkv(k).gt.winmax) then     !internal node
            inodv(k)=0
          else if(winkv(k).eq.0.) then  !out of system
            inodv(k)=-2
          else                          !boundary node
            inodv(k)=-1
          end if
        end do
c
c now mark open boundary nodes
c
	do ibc=1,nbc
          ibtyp = itybnd(ibc)
	  n = nkbnds(ibc)
          if(ibtyp.ge.3) then       !$$ibtyp3	!$$ibtyp4
            do ii=1,n
              k=kbnds(ibc,ii)
              if(inodv(k).eq.-1) inodv(k)=ibc
            end do
          else if(ibtyp.gt.0) then
            do ii=1,n
              k=kbnds(ibc,ii)
              if(inodv(k).ne.-1) goto 99
              inodv(k)=ibc
            end do
          end if
        end do
c
        return
   99   continue
        write(6,*) 'error for open boundary node'
        if( inodv(k) .eq. 0 ) then
          write(6,*) 'the node is an inner node'
        else if( inodv(k) .eq. -2 ) then
          write(6,*) 'the node is not in the system (dry)'
        else
          write(6,*) 'internal error (1)'
        end if
        write(6,*) 'k,inodv : ',ipext(k),inodv(k)
        stop 'error stop setnod : open boundary node'
        end
c
c****************************************************************

        function kthis(i,ie)

c gets node at position i in ie
c
c i     position
c ie    element

        implicit none

c arguments
	integer kthis
        integer i,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v

	kthis = nen3v(i,ie)

	end

c****************************************************************
c
        function knext(k,ie)
c
c gets node after k in ie
c
c k     actual node
c ie    element
c
        implicit none
c
c arguments
	integer knext
        integer k,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            knext=nen3v(mod(i,3)+1,ie)
            return
          end if
        end do
c
        knext=0
c
        return
        end
c
c****************************************************************
c
        function kbhnd(k,ie)
c
c gets node before k in ie
c
c k     actual node
c ie    element
c
        implicit none
c
c arguments
	integer kbhnd
        integer k,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            kbhnd=nen3v(mod(i+1,3)+1,ie)
            return
          end if
        end do
c
        kbhnd=0
c
        return
        end
c
c****************************************************************
c
        function ithis(k,ie)
c
c gets i of k in ie
c
c k     actual node
c ie    element
c
        implicit none
c
c arguments
	integer ithis
        integer k,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ithis=i
            return
          end if
        end do
c
        ithis=0
c
        return
        end
c
c****************************************************************
c
        function inext(k,ie)
c
c gets i after k in ie
c
c k     actual node
c ie    element
c
        implicit none
c
c arguments
	integer inext
        integer k,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            inext=mod(i,3)+1
            return
          end if
        end do
c
        inext=0
c
        return
        end
c
c****************************************************************
c
        function ibhnd(k,ie)
c
c gets i before k in ie
c
c k     actual node
c ie    element
c
        implicit none
c
c arguments
	integer ibhnd
        integer k,ie
c common
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer i
c
        do i=1,3
          if( nen3v(i,ie) .eq. k ) then
            ibhnd=mod(i+1,3)+1
            return
          end if
        end do
c
        ibhnd=0
c
        return
        end
c
c****************************************************************
c
        subroutine pntfla(k,ipf,ipl)
c
c gets first and last element in element link of node k
c
c absolute search without barene
c
c k     actual node
c ipf   first element (return)
c ipl   last  element (return)
c
        implicit none
c
c arguments
        integer k,ipf,ipl
c common
        integer ilinkv(1),lenkv(1)
        common /ilinkv/ilinkv,/lenkv/lenkv

        ipf=ilinkv(k)+1
        ipl=ilinkv(k+1)

	if( lenkv(ipl) .eq. 0 ) ipl = ipl - 1	!FIXME

        end
c
c****************************************************************
c
        subroutine pntfl(k,ipf,ipl)
c
c gets first and last element in element link of node k
c
c k     actual node
c ipf   first element (return)
c ipl   last  element (return)
c
        implicit none
c
c arguments
        integer k,ipf,ipl
c common
        integer ilinkv(1),lenkv(1),inodv(1)
        integer iwegv(1)
        common /ilinkv/ilinkv,/lenkv/lenkv,/inodv/inodv
        common /iwegv/iwegv
c local
        integer ie,ip,ip0,ip1
        logical bin
c functions
        logical iskbnd,iseins
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0
c
        if( iskbnd(k) ) then
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          ipf=0
          ipl=0
c
          bin=iseins(lenkv(ip1))
          do ip=ip0,ip1
            if( .not.bin .and. iseins(lenkv(ip)) ) then
              if(ipf.eq.0) then
                ipf=ip
              else  !here we have elements only attached by one vertex
                ipf=-1
                ipl=-1
                return
              end if
            end if
            bin=iseins(lenkv(ip))
          end do
c
          bin=iseins(lenkv(ip0))
          do ip=ip1,ip0,-1
            if( .not.bin .and. iseins(lenkv(ip)) ) then
              if(ipl.eq.0) then
                ipl=ip
              else
                ipf=-1
                ipl=-1
                return
              end if
            end if
            bin=iseins(lenkv(ip))
          end do
        else
          ipf=ilinkv(k)+1
          ipl=ilinkv(k+1)
        end if
c
        return
        end
c
c****************************************************************
c
        subroutine ntlink(ilinkv,lenkv)
c
c sets up vector with element links and a pointer to it
c
c ilinkv    pointer to links
c lenkv     link element numbers
c
c
c number of links of node n : nl = ilinkv(n+1)-ilinkv(n)
c links of node n           : ( lenkv ( ilinkv(n)+i ), i=1,nl )
c
        implicit none
c
c arguments
        integer ilinkv(1),lenkv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),inodv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /inodv/inodv
c local
        integer ie,i,n,k,k1,ip,ip0,ip1,ipe
c        integer nbn,nli
c        integer nli
c functions
        integer knext,kbhnd,ipext
	real winkk
        logical iskbnd
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
c
        call setweg(-1,n)
        call setnod
c
c first the total number of links for each node is established
c
        do k=1,nkn
          ilinkv(k)=0
        end do
c
        do ie=1,nel
          do i=1,3
            k=nen3v(i,ie)
            ilinkv(k)=ilinkv(k)+1
          end do
        end do
c
c boundary nodes have one link more --> find boundary nodes and add it
c
        do k=1,nkn
          if( iskbnd(k) ) ilinkv(k)=ilinkv(k)+1
        end do
c
c now create pointer into array
c
        n=0
        do k=1,nkn
          i=ilinkv(k)
          ilinkv(k)=n
          n=n+i
        end do
        ilinkv(nkn+1)=n
c
c now get the element numbers of the links
c
        do i=1,ilinkv(nkn+1)
          lenkv(i)=0
        end do
c
        do ie=1,nel
          do i=1,3
            k=nen3v(i,ie)
            ip=ilinkv(k)+1
            ip1=ilinkv(k+1)
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              ip=ip+1
            end do
            if(ip.gt.ip1) then  !error
              write(6,*) k,k1,ilinkv(k),ip,ip1
              stop 'error stop ntlink : not possible branch 1'
            end if
            lenkv(ip)=ie
          end do
        end do
c
c sort element entries
c
        do k=1,nkn
c
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
c
c if k is boundary node, find first element (in anti-clockwise sense)
c
          if( iskbnd(k) ) then
c
            ip=ip0
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              ipe=ip
              k1=knext(k,lenkv(ipe))
              ip=ip0
              do while(ip.le.ip1.and.lenkv(ip).ne.0)
                if(knext(k1,lenkv(ip)).eq.k) goto 1
                ip=ip+1
              end do
    1         continue
            end do
c
c           at ipe is first element --> now swap
c
            ie=lenkv(ip0)
            lenkv(ip0)=lenkv(ipe)
            lenkv(ipe)=ie
c
          end if
c
c sort next elements
c
          do while(ip0.le.ip1-1.and.lenkv(ip0).ne.0)
            k1=kbhnd(k,lenkv(ip0))
            ip=ip0+1
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              if(knext(k,lenkv(ip)).eq.k1) goto 2
              ip=ip+1
            end do
    2       continue
c
            ip0=ip0+1
c
            if(ip.le.ip1.and.lenkv(ip).ne.0) then  !swap here
              ie=lenkv(ip0)
              lenkv(ip0)=lenkv(ip)
              lenkv(ip)=ie
            end if
          end do
c
        end do  !loop over all nodes
c
c test internal and external nodes	!FIXME
c
        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)

	  do ip=ip0,ip1
	    if( lenkv(ip) .le. 0 ) then
c	      write(6,*) '0 link ',iskbnd(k),k,ip,ip0,ip1
	    end if
	  end do
	end do

        n=0
        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
c
          ip=ip0+1
          do while(ip.le.ip1.and.lenkv(ip).ne.0)
            if(kbhnd(k,lenkv(ip-1)).ne.knext(k,lenkv(ip))) then
              n=n+1
              write(6,*) 'node ',ipext(k),' has angle ',winkk(k)
            end if
            ip=ip+1
          end do
        end do
c
        if( n.ne.0 ) then
          write(6,*) 'Nodes marked as internal are boundary nodes.'
          write(6,*) 'Be sure to set the parameter WINMAX'
          write(6,*) 'in routine SETNOD (file SUBLNK)'
          write(6,*) 'higher than maximum angle of nodes above'
          write(6,*) 'but lower than 360.'
          stop 'error stop ntlink'
        end if
c
        return
        end
c
c****************************************************************
c
        subroutine setkan(kantv)
c
c sets up vector kantv
c
        implicit none
c
c arguments
        integer kantv(2,1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1), inodv(1)
        integer iwegv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv, /inodv/inodv
        common /iwegv/iwegv
c local
c       integer k,ie,ipf,ipl !$$ALPHA
        integer k,ipf,ipl
c functions
        integer knext,kbhnd
        logical iskbnd
c	logical iseins	!$$ALPHA
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
c       iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0	!$$ALPHA
c
        do k=1,nkn
          if( iskbnd(k) ) then
            call pntfl(k,ipf,ipl)
            if(ipf.eq.-1) then
              kantv(1,k)=0
              kantv(2,k)=0
            else
              kantv(1,k) = knext( k , lenkv(ipf) )
              kantv(2,k) = kbhnd( k , lenkv(ipl) )
            end if
          else
            kantv(1,k)=0
            kantv(2,k)=0
          end if
        end do
c
        return
        end
c
c****************************************************************
c
        subroutine setdxy(dxv,dyv)
c
c sets up vector dxv and dyv
c
        implicit none
c
c arguments
        real dxv(1),dyv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1), inodv(1)
        integer iwegv(1)
        real xgv(1),ygv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv, /inodv/inodv
        common /iwegv/iwegv
        common /xgv/xgv, /ygv/ygv
c local
c       integer k,k1,k2,ie,ipf,ipl !$$ALPHA
        integer k,k1,k2,ipf,ipl
c functions
        integer knext,kbhnd
        logical iskbnd,isobnd
c       logical iseins	!$$ALPHA
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
        isobnd(k) = inodv(k).gt.0
c       iseins(ie) = ie.gt.0.and.iwegv(ie).eq.0	!$$ALPHA
c
        do k=1,nkn
          if( iskbnd(k) ) then                    !boundary node
            call pntfl(k,ipf,ipl)
            if(ipf.eq.-1) then                      !irregular boundary
              dxv(k)=0.
              dyv(k)=0.
            else                                    !regular
              k1 = knext( k , lenkv(ipf) )
              k2 = kbhnd( k , lenkv(ipl) )
              if( isobnd(k) ) then                    !open boundary node
                if( .not. isobnd(k1) ) then             !k1 is not open
                  dxv(k) = xgv(k1)-xgv(k)
                  dyv(k) = ygv(k1)-ygv(k)
                else if ( .not. isobnd(k2) ) then       !k2 is not open
                  dxv(k) = xgv(k2)-xgv(k)
                  dyv(k) = ygv(k2)-ygv(k)
                else                                    !central open
                  dyv(k) =    xgv(k1)-xgv(k2)
                  dxv(k) = -( ygv(k1)-ygv(k2) )
                end if
              else                                    !normal boundary node
                dxv(k)=xgv(k1)-xgv(k2)
                dyv(k)=ygv(k1)-ygv(k2)
              end if
            end if
          else                                    !no boundary node
            dxv(k)=0.
            dyv(k)=0.
          end if
        end do
c
        return
        end

c
c****************************************************************
c****************************************************************
c****************************************************************
c
        subroutine setelt
c
c sets up vector ieltv
c
        implicit none
c
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1), inodv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv, /inodv/inodv
        integer ieltv(3,1)
        common /ieltv/ieltv
c local
        integer k,ie,ii,ip,ip0,ip1
c functions
        integer inext,kbhnd
        logical isobnd
        isobnd(k) = inodv(k).gt.0
c
        do ie=1,nel
          do ii=1,3
            ieltv(ii,ie)=0
          end do
        end do
c
        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)-1
c
          do ip=ip0,ip1
            ie=lenkv(ip)
            ieltv(inext(k,ie),ie)=lenkv(ip+1)
          end do
c
c last element
c
          ie=lenkv(ip)
          if(ie.gt.0) then                              !is internal node
            ieltv(inext(k,ie),ie)=lenkv(ip0)
          else                                          !is boundary node
            ie=lenkv(ip-1)
            if( isobnd(k) .and. isobnd(kbhnd(k,ie)) ) then  !open boundary
              ieltv(inext(k,ie),ie)=-1
            end if
          end if
        end do
c
        return
        end
c
c****************************************************************

	subroutine checkelt

c checks integrity of ieltv

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1), inodv(1)
        common /ilinkv/ilinkv, /lenkv/lenkv, /inodv/inodv
        integer ieltv(3,1)
        common /ieltv/ieltv
        integer nen3v(3,1)
        common /nen3v/nen3v
c local
        integer ie,ii,ien,ienn,inn,k,kn
	integer nbnd,nobnd,nintern
c functions
        integer inext,knext

	nbnd = 0
	nobnd = 0
	nintern = 0

        do ie=1,nel
          do ii=1,3
	    nintern = nintern + 1
            ien = ieltv(ii,ie)
	    if( ien .gt. 0 ) then
		if( ien .gt. nel ) goto 99
	    	k = nen3v(ii,ie)
	    	kn = knext(k,ie)
		inn = inext(kn,ien)
		ienn = ieltv(inn,ien)
		if( ie .ne. ienn ) goto 98
	    else if( ien .eq. 0 ) then
		nbnd = nbnd + 1
	    else if( ien .eq. -1 ) then
		nobnd = nobnd + 1
	    else
		goto 99
	    end if
          end do
	end do

	write(6,*) 'checkelt is ok'
	write(6,*) '  internal sides =      ',nintern
	write(6,*) '  boundary sides =      ',nbnd
	write(6,*) '  open boundary sides = ',nobnd
	write(6,*) '  total sides =         ',nintern+nbnd+nobnd

	return
   98	continue
	write(6,*) 'ie,ii,ien,nel: ',ie,ii,ien,nel
	write(6,*) 'k,kn,inn,ienn: ',k,kn,inn,ienn
	write(6,*) 'nen3v: ',ie,(nen3v(ii,ie),ii=1,3)
	write(6,*) 'nen3v: ',ien,(nen3v(ii,ien),ii=1,3)
	write(6,*) 'ieltv: ',ie,(ieltv(ii,ie),ii=1,3)
	write(6,*) 'ieltv: ',ien,(ieltv(ii,ien),ii=1,3)
	stop 'error stop checkelt: corrupt data structure of ieltv'
   99	continue
	write(6,*) 'ie,ii,ien,nel: ',ie,ii,ien,nel
	stop 'error stop checkelt: corrupt data structure of ieltv'
	end

c****************************************************************

        function ielt(i,ie)

c gets i'th neigbour of element ie, i in [1-3]

        implicit none

	integer ielt
        integer i,ie

        integer ieltv(3,1)
        common /ieltv/ieltv

	ielt = ieltv(i,ie)

        end

c****************************************************************
c****************************************************************
c****************************************************************
c
        subroutine setlin(linkv)
c
c sets up vector linkv
c
        implicit none
c
c arguments
        integer linkv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv
c local
        integer k,ip,ip0,ip1
c functions
        integer knext,kbhnd
c
        do k=1,nkn
c
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
c
          do ip=ip0,ip1-1
            linkv(ip)=knext(k,lenkv(ip))
          end do
c
          if(lenkv(ip1).gt.0) then  !no boundary node
            linkv(ip1)=knext(k,lenkv(ip1))
          else
            linkv(ip1)=kbhnd(k,lenkv(ip1-1))
          end if
c
        end do
c
        return
        end
c
c****************************************************************
c
        subroutine getlin(k,link,nlink)
c
c gets links for node k
c
c k       node to links
c link    links to node k (return)
c nlink   total number of links to node k (return)
c
c link must be big enough (ngr)
c
        implicit none
c
c arguments
        integer k,link(1),nlink
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv
c local
        integer ip,ip0,ip1
c functions
        integer knext,kbhnd
c
        ip0=ilinkv(k)+1
        ip1=ilinkv(k+1)
        nlink=0
c
        do ip=ip0,ip1-1
          nlink=nlink+1
          link(nlink)=knext(k,lenkv(ip))
        end do
c
        nlink=nlink+1
        if(lenkv(ip1).gt.0) then  !no boundary node
          link(nlink)=knext(k,lenkv(ip1))
        else
          link(nlink)=kbhnd(k,lenkv(ip1-1))
        end if
c
        return
        end
c
c****************************************************************
c
        subroutine tstlnk
c
c sets up vector with element links and a pointer to it
c
        implicit none
c
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),inodv(1)
        integer ilinkv(1),lenkv(1),linkv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /inodv/inodv
        common /ilinkv/ilinkv, /lenkv/lenkv, /linkv/linkv
        integer ilink1v(1),lenk1v(1),link1v(1)
        common /ilink1v/ilink1v, /lenk1v/lenk1v, /link1v/link1v
        integer kantv(2,1)
        integer kant1v(2,1)
        common /kantv/kantv
        common /kant1v/kant1v
	integer ieltv(3,1)
	common /ieltv/ieltv
	integer ielt1v(3,1)
	common /ielt1v/ielt1v
c local
	logical bstop
        integer ie,i,j,n,k,ip,ip0,ip1,ii
        integer nli
	integer nlidim
c       integer nnkn,nnel,nnbn,nnli,nnis,nnod
c functions
        integer knext,kbhnd
        logical iskbnd
        iskbnd(k) = inodv(k).ne.0 .and. inodv(k).ne.-2
c
        call ntlink(ilinkv,lenkv)
	call getdim('nlidim',nlidim)
	nlidim = 2 * nlidim
	write(6,*) 'setting up link structure... ',nlidim
	call mklenk(nlidim,nkn,nel,nen3v,ilink1v,lenk1v)

	call setlin(linkv)
	call mklink(nkn,ilink1v,lenk1v,link1v)

	write(6,*) 'checking ilink structure... '
	bstop = .false.
	do k=1,nkn+1
	  if( ilinkv(k) .ne. ilink1v(k) ) then
		write(6,*) k,ilinkv(k),ilink1v(k)
		bstop = .true.
	  end if
	end do
	write(6,*) 'checking lenk structure... ',ilinkv(nkn+1)
	do k=1,ilinkv(nkn+1)
	  if( lenkv(k) .ne. lenk1v(k) ) then
		write(6,*) k,lenkv(k),lenk1v(k)
		bstop = .true.
	  end if
	end do
	write(6,*) 'checking link structure... ',ilinkv(nkn+1)
	do k=1,ilinkv(nkn+1)
	  if( linkv(k) .ne. link1v(k) ) then
		write(6,*) k,linkv(k),link1v(k)
		bstop = .true.
	  end if
	end do
	if( bstop ) stop 'error stop tstlnk'

	call setkan(kantv)
	call mkkant(nkn,ilinkv,lenkv,linkv,kant1v)
	write(6,*) '--- 1'
	call checkkant(nkn,kantv)
	write(6,*) '--- 2'
	call checkkant(nkn,kant1v)
	write(6,*) '--- 3'
	do k=1,nkn
	  if( kantv(1,k) .ne. kant1v(1,k) ) then
		write(6,*) k,1,kantv(1,k),kant1v(1,k)
		bstop = .true.
	  end if
	  if( kantv(2,k) .ne. kant1v(2,k) ) then
		write(6,*) k,2,kantv(2,k),kant1v(2,k)
		bstop = .true.
	  end if
	end do
	if( bstop ) stop 'error stop tstlnk (2)'

        call setweg(-1,n)
        call setnod
c
c we also set ieltv (may be not needed)
c
        call setelt
	call checkelt
	call mkielt(nkn,nel,ilinkv,lenkv,linkv,ielt1v)
	call checkielt(nel,ieltv)
	call checkielt(nel,ielt1v)
	write(6,*) 'checkielt ok...'
	call update_ielt(nel,inodv,ielt1v)
	write(6,*) 'ieltv updated...'
	do ie=1,nel
	 do ii=1,3
	  if( ieltv(ii,ie) .ne. ielt1v(ii,ie) ) then
		write(6,*) ie,ii,ieltv(ii,ie),ielt1v(ii,ie)
		bstop = .true.
	  end if
	 end do
	end do
	if( bstop ) stop 'error stop tstlnk (3)'
c
c write some statistics
c
        nli=ilinkv(nkn+1)/2
        n=0
        do i=1,2*nli
          if(lenkv(i).eq.0) n=n+1
        end do
c
        write(6,*) 'nel                      : ',nel
        write(6,*) 'links from ilink(pointer): ',ilinkv(nkn+1)
        write(6,*) 'double links             : ',2*nli
        write(6,*) 'links                    : ',nli
        write(6,*) 'links not set            : ',n
c
c control array link
c
        n=0
        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          ip=ip0-1
          do while(ip.lt.ip1)
            ip=ip+1
            ie=lenkv(ip)
            i=ip+1
            do while(i.lt.ip1)
              if(lenkv(i).eq.ie) then
                write(6,*) 'elements are not unique'
                write(6,*) 'k,ip,ip1 : ',k,ip,ip1
                write(6,*) (lenkv(j),j=ip0,ip1)
                stop 'error stop : tstlnk'
              end if
              n=n+1
              i=i+1
            end do
          end do
        end do
c
c        write(6,*) 'comparisons made (unique): ',n
c
        n=0
        do i=1,2*nli
          if(lenkv(i).eq.0) n=n+1
        end do
c
c        write(6,*) 'number of nodes not set  : ',n
c
        n=0
        do k=1,nkn
c
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
c
c test sorting
c
          ip=ip0+1
          do while(ip.le.ip1.and.lenkv(ip).ne.0)
            if(kbhnd(k,lenkv(ip-1)).ne.knext(k,lenkv(ip))) then
              write(6,*) 'elements not sorted'
              write(6,*) 'k,ip : ',k,ip
              write(6,*) (lenkv(j),j=ip0,ip1)
              stop 'error stop : tstlnk'
            end if
            n=n+1
            ip=ip+1
          end do
c
c test last element if 0
c
          if(lenkv(ip1).eq.0.and..not.iskbnd(k)) then
            write(6,*) 'element 0 and no boundary node'
            write(6,*) 'k : ',k
            write(6,*) (lenkv(j),j=ip0,ip1)
            stop 'error stop : tstlnk'
          end if
c
c test for more elements 0
c
          if(lenkv(ip1).eq.0) then
            do ip=ip0,ip1-1
              if(lenkv(ip).eq.0) then
                write(6,*) 'more than one element 0'
                write(6,*) 'k : ',k
                write(6,*) (lenkv(j),j=ip0,ip1)
                stop 'error stop : tstlnk'
              end if
            end do
          end if
c
        end do
c
c        write(6,*) 'comparisons made (sort)  : ',n
c
c        call nknel(nnkn,nnel,nnbn,nnli,nnis,nnod)	!$$NKNEL
c
c        write(6,*) 'nkn,nel : ',nnkn,nnel
c        write(6,*) 'nbn,nli : ',nnbn,nnli
c        write(6,*) 'nis,nod : ',nnis,nnod
c
        return
        end
c
c*****************************************************************

        subroutine newlnk

c administrates new topological routines
c
c ... iwegv has already been set

        implicit none

c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        integer ieltv(3,1),kantv(2,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
        common /ieltv/ieltv,/kantv/kantv
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

        call setnar(nnar)
        call setwnk(wink)

        nnnis = iround( (wink-nnbn*180.)/360. ) + nnar

        nnnel=2*nnkn-nnbn+2*(nnnis-nnar)-nnod

        write(n88,'(a,i8,10i6)') 'newlnk: ',it
     +          ,nnkn,nnel,nnbn,nnli,nnis,nnod,nnar
     +          ,nnnis,nnnel-nnel,nel-nnel

        end

c****************************************************************
c
        function winkk(k)
c
c total angle at node k
c
c k       node to links
c
        implicit none
c
c arguments
	real winkk
        integer k
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilinkv(1), lenkv(1), iwegv(1)
        real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /ilinkv/ilinkv, /lenkv/lenkv, /iwegv/iwegv
        common /ev/ev
c local
        integer ip,ip0,ip1,ie
        real w
c functions
        integer ithis
c
        ip0=ilinkv(k)+1
        ip1=ilinkv(k+1)
c
        w=0.
        do ip=ip0,ip1
          ie=lenkv(ip)
          if( ie.gt.0 .and. iwegv(ie).eq.0 ) then
            w = w + ev(10+ithis(k,ie),ie)
          end if
        end do
c
        winkk=w
c
        return
        end
c
c****************************************************************
c
        subroutine setnar(nar)
c
c sets number of areas
c
c nar   number of areas (return)
c
        implicit none
c
c arguments
        integer nar
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),ieltv(3,1),iwegv(1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /ieltv/ieltv, /iwegv/iwegv
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
              if(ien.gt.0.and.iwegv(ien).lt.3) then !possible to enter
                if(iwegv(ien).eq.0.or.iwegv(ieo).eq.3) then !may enter
                  i=inext( nen3v(mod(i,3)+1,ieo) , ien )
                  ieo=ien
                  n2=n2+1
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
        implicit none
c
c arguments
        real wink
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1),inodv(1),iwegv(1)
        real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nen3v/nen3v, /inodv/inodv, /iwegv/iwegv
        common /ev/ev
c local
        integer ie,ii,k
        real w
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
        implicit none
c
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        integer iwegv(1)
        real ev(13,1)
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /femtim/ itanf,itend,idt,nits,niter,it
        common /iwegv/iwegv
        common /ev/ev
c local
        real arin,arout,artot,area,rw
        integer ie
c functions
	real zvbnds
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

        if(nbc.gt.0) then
          rw=zvbnds(1)
        else
          rw=0.
        end if

c        write(88,'(i8,f12.3,2f12.2,e12.4)') it,rw,arin,arout,artot
c
        return
        end

c****************************************************************
c****************************************************************
c****************************************************************

	function is_internal_node(k)

	implicit none

	logical is_internal_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_internal_node = inodv(k) .eq. 0

	end

c****************************************************************

	function is_boundary_node(k)

	implicit none

	logical is_boundary_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_boundary_node = inodv(k) .ne. 0 .and. inodv(k) .ne. -2

	end

c****************************************************************

	function is_open_boundary_node(k)

	implicit none

	logical is_open_boundary_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_open_boundary_node = inodv(k) .gt. 0

	end

c****************************************************************

	function is_dry_node(k)

	implicit none

	logical is_dry_node
	integer k

        integer inodv(1)
        common /inodv/inodv

	is_dry_node = inodv(k) .eq. -2

	end

c****************************************************************
c****************************************************************
c****************************************************************

