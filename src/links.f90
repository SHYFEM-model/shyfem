!
! $Id: links.f,v 1.8 2009-05-21 09:24:00 georg Exp $
!
! geometric routines (static)
!
! contents :
!
! revision log :
!
! 10.08.2003    ggu	new routines for ieltv, kantv
! 18.10.2005    ggu	more debug output, new routines get_link, get_lenk
! 06.04.2009    ggu	changed nlidim to nlkddi
! 30.03.2011    ggu	bverbose introduced
! 30.05.2015    ggu	new subroutine mklenkii to compute lenkiiv
! 02.12.2015    ggu	routines get_link() and get_lenk() deleted
!
!****************************************************************
!****************************************************************
!****************************************************************

      module links

        use lnku
        contains

!****************************************************************
        subroutine mklenk(nlkddi,nkn,nel,nen3v,ilinkv,lenkv)

! sets up vector with element links and a pointer to it
!
! only array nen3v is needed, no aux array is needed
!
! ilinkv    pointer to links
! lenkv     link to element numbers
!
! number of links of node n : nl = ilinkv(n+1)-ilinkv(n)
! links of node n           : ( lenkv ( ilinkv(n)+i ), i=1,nl )

        implicit none

! arguments
        integer nlkddi
        integer nkn
        integer nel
        integer nen3v(3,nel)
        integer ilinkv(nkn+1)
        integer lenkv(*)
! local
        logical binside
        integer nloop
        integer ie,i,n,k
        integer ip,ip0,ip1,ipe
        integer ie0,ie1,k0,k1
        integer nfill
! functions

	ip = 0
	ip1 = 0

!-------------------------------------------------------------------
! first the total number of elements (links) for each node is established
!-------------------------------------------------------------------

        do k=1,nkn
          ilinkv(k)=1           !one more link -> will be corrected later
        end do

        do ie=1,nel
          do i=1,3
            k=nen3v(i,ie)
            ilinkv(k)=ilinkv(k)+1
          end do
        end do

!-------------------------------------------------------------------
! now create pointer into array
!-------------------------------------------------------------------

        n=0
        do k=1,nkn
          i=ilinkv(k)
          ilinkv(k)=n
          n=n+i
        end do
        ilinkv(nkn+1)=n

        if( n .gt. nlkddi ) goto 98

!-------------------------------------------------------------------
! now enter the element numbers of the links
!-------------------------------------------------------------------

        do i=1,ilinkv(nkn+1)
          lenkv(i)=0
        end do

        do ie=1,nel
          do i=1,3
            k=nen3v(i,ie)
	    if( k .gt. nkn ) goto 97
            ip=ilinkv(k)+1      !first entry
            ip1=ilinkv(k+1)     !last possible entry
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              ip=ip+1
            end do
            if(ip.gt.ip1) goto 99       !error 1
            lenkv(ip)=ie
          end do
        end do

!-------------------------------------------------------------------
! sort element entries
!-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)

!         ----------------------------------------------------------
!         if k is boundary node, find first element (in anti-clockwise sense)
!         ----------------------------------------------------------

          binside = .true.
          nloop = 0
          ip=ip0

          do while( binside .and. nloop .le. ip1-ip0 )
            ipe=ip
            k1=knext(k,lenkv(ipe))
            ip=ip0
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              if(knext(k1,lenkv(ip)).eq.k) goto 1
              ip=ip+1
            end do
            binside = .false.
    1       continue
            nloop = nloop + 1
          end do

!         --------------------------------------------------------
!         at ipe is first element --> now swap
!         --------------------------------------------------------

          if( .not. binside ) then      !boundary element
            ie=lenkv(ip0)
            lenkv(ip0)=lenkv(ipe)
            lenkv(ipe)=ie
          end if

!         ----------------------------------------------------------
!         sort next elements
!         ----------------------------------------------------------

          do while(ip0.le.ip1-1.and.lenkv(ip0).ne.0)
            k1=kbhnd(k,lenkv(ip0))
            ip=ip0+1
            do while(ip.le.ip1.and.lenkv(ip).ne.0)
              if(knext(k,lenkv(ip)).eq.k1) goto 2
              ip=ip+1
            end do
    2       continue

            ip0=ip0+1

            if(ip.le.ip1.and.lenkv(ip).ne.0) then  !swap here
              ie=lenkv(ip0)
              lenkv(ip0)=lenkv(ip)
              lenkv(ip)=ie
            end if
          end do

        end do  !loop over all nodes

!-------------------------------------------------------------------
! compact structure of non boundary nodes
!
! for elements we could theoretically compact all entries,
! but we want to reuse ilinkv also for the node pointer,
! and therefore we need one more entry for boundary nodes
!-------------------------------------------------------------------

        nfill = 0

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

          ie0 = lenkv(ip0)
          ie1 = lenkv(ip1)
          k0 = knext(k,ie0)
          k1 = kbhnd(k,ie1)
          if( k0 .ne. k1 ) ip1 = ip1 + 1        !boundary node

          ilinkv(k) = nfill
          do ip=ip0,ip1
            nfill = nfill + 1
            lenkv(nfill) = lenkv(ip)
          end do

        end do

        ilinkv(nkn+1) = nfill

!-------------------------------------------------------------------
! check structure
!-------------------------------------------------------------------

	call checklenk(nkn,ilinkv,lenkv,nen3v)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        return
   97   continue
        write(6,*) 'node: ',k,'  nkn: ',nkn
        write(6,*) 'error stop mklenk : internal error (2)'
        stop 'error stop mklenk : internal error (2)'
   98   continue
        write(6,*) n,nlkddi
        write(6,*) nkn,nel,2*nkn+nel
        stop 'error stop mklenk: nlkddi'
   99   continue
        write(6,*) 'node: ',k
        write(6,*) 'element: ',ie
        write(6,*) 'first entry: ',ilinkv(k)+1
        write(6,*) 'last entry: ',ip1
        write(6,*) 'actual pointer: ',ip
        write(6,*)'error stop mklenk : internal error (1)'
        stop 'error stop mklenk : internal error (1)'
        end

!****************************************************************

        subroutine checklenk(nkn,ilinkv,lenkv,nen3v)

! checks structure of lenkv
        implicit none

! arguments
        integer nkn
        integer ilinkv(nkn+1)
        integer lenkv(*)
	integer nen3v(3,*)
! local
	integer nbnd,nnull
	integer k,k0,k1,ie,ie0,ie1
	integer ip,ip0,ip1
	integer i
	integer ipp,ii
	logical bverbose

	bverbose = .false.

        nbnd = 0        !total number of boundary nodes
        nnull = 0       !total number of 0 entries

!-------------------------------------------------------------------
! loop over nodes
!-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

          do ip=ip0,ip1
            ie = lenkv(ip)
            if( ie .eq. 0 ) then	!0 in index found
              write(6,*) 'Node (internal) k = ',k
              write(6,*) k,ip0,ip1
              write(6,*) (lenkv(i),i=ip0,ip1)
              write(6,*)'error stop checklenk: structure of lenkv(2)'
              stop 'error stop checklenk: structure of lenkv (2)'
            end if
          end do

          do ip=ip0+1,ip1
            ie0 = lenkv(ip-1)
            ie1 = lenkv(ip)
            k0 = kbhnd(k,ie0)
            k1 = knext(k,ie1)
            if( k0 .ne. k1 ) then	!something is wrong
              write(6,*) 'Node (internal) k = ',k
              write(6,*) ip0,ip1,ip
              write(6,*) ie0,ie1,k0,k1
              write(6,*) ie0,(kthis(i,ie0),i=1,3)
              write(6,*) ie1,(kthis(i,ie1),i=1,3)
              write(6,*) (lenkv(i),i=ip0,ip1)
	      do ipp=ip0,ip1
		ie = lenkv(ipp)
	        write(6,*) ie,(nen3v(ii,ie),ii=1,3)
	      end do
              write(6,*)'error stop checklenk: structure of lenkv(3)'
              stop 'error stop checklenk: structure of lenkv(3)'
            end if
          end do

          ie0 = lenkv(ip0)
          ie1 = lenkv(ip1)
          k0 = knext(k,ie0)
          k1 = kbhnd(k,ie1)
          if( k0 .ne. k1 ) nbnd = nbnd + 1

          if( lenkv(ilinkv(k+1)) .eq. 0 ) nnull = nnull + 1

        end do

	if( bverbose ) then
          write(6,*) 'checklenk: ',nnull,nbnd,nkn
	end if

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine mklenkii(nlkddi,nkn,nel,nen3v,ilinkv,lenkv,lenkiiv)

! sets up vector lenkiiv
!
! only array nen3v is needed, no aux array is needed
!
! ilinkv    pointer to links
! lenkv     link to element numbers
!
! number of links of node n : nl = ilinkv(n+1)-ilinkv(n)
! links of node n           : ( lenkv ( ilinkv(n)+i ), i=1,nl )

        implicit none

! arguments
        integer nlkddi
        integer nkn
        integer nel
        integer nen3v(3,nel)
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer lenkiiv(*)
! local
        integer ie,i,n,k,ii,ibase

!-------------------------------------------------------------------
! compute the vertec number for elements connected to node k
!-------------------------------------------------------------------

	do i=1,nlkddi
	  lenkiiv(i) = 0
	end do

        do k=1,nkn
	  n = ilinkv(k+1)-ilinkv(k)
	  ibase = ilinkv(k)
	  if( lenkv(ibase+n) .le. 0 ) n = n - 1
	  do i=1,n
	    ie = lenkv(ibase+i)
	    do ii=1,3
	      if( nen3v(ii,ie) == k ) exit
	    end do
	    if( ii > 3 ) goto 99
	    lenkiiv(ibase+i) = ii
	  end do
        end do

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	return
   99	continue
	write(6,*) 'k,ie,n,ibase ',k,ie,n,ibase
	write(6,*)'error stop mklenkii: internal error (1)'
	stop 'error stop mklenkii: internal error (1)'
	end

!****************************************************************

        subroutine mklink(nkn,ilinkv,lenkv,linkv)

! sets up vector with node links
!
! ilinkv and lenkv must have already been set up
!
! only array nen3v is needed, no aux array is needed
!
! ilinkv    pointer to links
! lenkv     link to element numbers
! linkv     link to node numbers

        implicit none

! arguments
        integer nkn
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
! local
        integer ie,k
        integer ip,ip0,ip1
! functions

!-------------------------------------------------------------------
! loop over nodes
!-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) ip1 = ip1 - 1

          do ip=ip0,ip1
            ie = lenkv(ip)
            linkv(ip) = knext(k,ie)
          end do

!         ----------------------------------------------------------
!         if boundary node insert last node
!         ----------------------------------------------------------

          ip=ilinkv(k+1)
          if( lenkv(ip) .eq. 0 ) then
            linkv(ip) = kbhnd(k,ie)
          end if

        end do

!-------------------------------------------------------------------
! check structure
!-------------------------------------------------------------------

	call checklink(nkn,ilinkv,linkv)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        end

!****************************************************************

        subroutine checklink(nkn,ilinkv,linkv)

! checks structure of lenkv

        implicit none

! arguments
        integer nkn
        integer ilinkv(nkn+1)
        integer linkv(*)
! local
	integer k,k1,i
	integer ip,ip0,ip1
	integer ipk,ipk0,ipk1
	logical bverbose

	bverbose = .false.

	k1 = 0				!compiler warnings
	ipk0 = 0
	ipk1 = 0

!-------------------------------------------------------------------
! loop over nodes
!-------------------------------------------------------------------

        do k=1,nkn

          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( linkv(ip1) .eq. 0 ) then	!0 in index found
            write(6,*) 'Node (internal) k = ',k
            write(6,*) ip0,ip1
            write(6,*) (linkv(ip),ip=ip0,ip1)
            write(6,*)'error stop checklink: internal error (1)'
            stop 'error stop checklink: internal error (1)'
          end if

	  do ip=ip0,ip1
	    k1 = linkv(ip)
            ipk0=ilinkv(k1)+1
            ipk1=ilinkv(k1+1)
	    ipk = ipk0
	    do while( ipk .le. ipk1 .and. linkv(ipk) .ne. k )
	      ipk = ipk + 1
	    end do
	    if( ipk .gt. ipk1 ) then	!node not found
              write(6,*) 'Node (internal) k = ',k
              write(6,*) ip0,ip1
              write(6,*) (linkv(i),i=ip0,ip1)
              write(6,*) k1,ipk0,ipk1
              write(6,*) (linkv(i),i=ipk0,ipk1)
              write(6,*)'error stop checklink: internal error (2)'
              stop 'error stop checklink: internal error (2)'
	    end if
	  end do
        end do

	if( bverbose ) then
          write(6,*) 'checklink: ',nkn
	end if

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        end

!****************************************************************

        subroutine mkkant(nkn,ilinkv,lenkv,linkv,kantv)

! makes vector kantv

        implicit none

! arguments
        integer nkn
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
        integer kantv(2,nkn)
! local
	integer k
	integer ip,ip0,ip1

        do k=1,nkn
          kantv(1,k)=0
          kantv(2,k)=0
        end do

        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1)
          if( lenkv(ip1) .eq. 0 ) then		!boundary node
            kantv(1,k) = linkv(ip0)
            kantv(2,k) = linkv(ip1)
          end if
        end do

!-------------------------------------------------------------------
! check structure
!-------------------------------------------------------------------

	call checkkant(nkn,kantv)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        end

!****************************************************************

        subroutine checkkant(nkn,kantv)

! checks structure of kantv

        implicit none

! arguments
        integer nkn
        integer kantv(2,nkn)
! local
	integer k,k1,k2
	integer nbnd,nint
	logical bverbose

	bverbose = .false.

	nbnd = 0
	nint = 0

!-------------------------------------------------------------------
! loop over nodes
!-------------------------------------------------------------------

        do k=1,nkn
	  k1 = kantv(1,k)
	  k2 = kantv(2,k)
	  if( k1 .gt. 0 .and. k2 .gt. 0 ) then
	    nbnd = nbnd + 1
	    if( k .ne. kantv(2,k1) .or. k .ne. kantv(1,k2) ) then
              write(6,*) 'Node (internal) k = ',k
	      write(6,*) 'k1,k2: ',k1,k2
	      write(6,*) 'backlink: ',kantv(2,k1),kantv(1,k2)
	      write(6,*) 'error stop checkkant: structure of kantv(2)'
	      stop 'error stop checkkant: structure of kantv(2)'
	    end if
	  else if( k1 .eq. 0 .and. k2 .eq. 0 ) then
	    nint = nint + 1
	  else
            write(6,*) 'Node (internal) k = ',k
	    write(6,*) 'k1,k2: ',k1,k2
	    write(6,*) 'error stop checkkant: structure of kantv(1)'
	    stop 'error stop checkkant: structure of kantv(1)'
	  end if
        end do

	if( bverbose ) then
	  write(6,*) 'checkkant: ',nkn,nint,nbnd
	end if

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        end

!****************************************************************

        subroutine mkielt(nkn,nel,ilinkv,lenkv,linkv,ieltv)

! makes vector ieltv (without open boundary nodes)

        implicit none

! arguments
        integer nkn,nel
        integer ilinkv(nkn+1)
        integer lenkv(*)
        integer linkv(*)
        integer ieltv(3,nel)
! local
	integer k,ie,ii
	integer ip,ip0,ip1

!-------------------------------------------------------------------
! initialize
!-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
            ieltv(ii,ie) = 0
          end do
        end do

!-------------------------------------------------------------------
! loop over links
!-------------------------------------------------------------------

        do k=1,nkn
          ip0=ilinkv(k)+1
          ip1=ilinkv(k+1) - 1

          do ip=ip0,ip1
            ie=lenkv(ip)
            ieltv(inext(k,ie),ie)=lenkv(ip+1)
          end do

!	  -------------------------------------------------------------------
!	  last element
!	  -------------------------------------------------------------------

          ie=lenkv(ip1+1)
          if(ie.gt.0) then                              !is internal node
            ieltv(inext(k,ie),ie)=lenkv(ip0)
          end if

        end do

!-------------------------------------------------------------------
! check structure
!-------------------------------------------------------------------
 
	call checkielt(nel,ieltv)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	end

!****************************************************************

        subroutine checkielt(nel,ieltv)

! checks structure of ieltv

        implicit none

! arguments
        integer nel
        integer ieltv(3,nel)
! local
	integer k,ie,ii
	integer kn,inn,ien,ienn
        integer nbnd,nobnd,nintern
	logical bverbose

!-------------------------------------------------------------------
! initialize
!-------------------------------------------------------------------

	bverbose = .false.

        nbnd = 0
        nobnd = 0
        nintern = 0

!-------------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
            ien = ieltv(ii,ie)
            if( ien .gt. 0 ) then
            	nintern = nintern + 1
                if( ien .gt. nel ) goto 99
                k = kthis(ii,ie)
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

	if( bverbose ) then
          write(6,*) 'checkielt is ok'
          write(6,*) '  internal sides =      ',nintern
          write(6,*) '  boundary sides =      ',nbnd
          write(6,*) '  open boundary sides = ',nobnd
          write(6,*) '  total sides =         ',nintern+nbnd+nobnd
	end if

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

        return
   98   continue
          write(6,*) 'Element (internal) ie = ',ie
          write(6,*) 'ii,ien,nel: ',ii,ien,nel
          write(6,*) 'k,kn,inn,ienn: ',k,kn,inn,ienn
          write(6,*) 'nen3v: ',ie,(kthis(ii,ie),ii=1,3)
          write(6,*) 'nen3v: ',ien,(kthis(ii,ien),ii=1,3)
          write(6,*) 'ieltv: ',ie,(ieltv(ii,ie),ii=1,3)
          write(6,*) 'ieltv: ',ien,(ieltv(ii,ien),ii=1,3)
          write(6,*)'error checkielt:corrupt data structure of ieltv(1)'
          stop 'error stop checkielt:corrupt data structure of ieltv(1)'
   99   continue
          write(6,*) 'Element (internal) ie = ',ie
          write(6,*) 'ii,ien,nel: ',ii,ien,nel
          write(6,*)'error checkielt:corrupt data structure of ieltv(2)'
          stop 'error stop checkielt:corrupt data structure of ieltv(2)'
	end

!****************************************************************

        subroutine update_ielt(nel,ibound,ieltv)

! updates vector ieltv with open boundary nodes

        implicit none

! arguments
        integer nel
        integer ibound(*)	! >0 => open boundary node
        integer ieltv(3,nel)
! local
	integer k,ie,ii,i
	integer ip,ip0,ip1

!-------------------------------------------------------------------
! loop over elements
!-------------------------------------------------------------------

        do ie=1,nel
          do ii=1,3
	    k = kthis(ii,ie)
	    if( ibound(k) .gt. 0 .and. ibound(knext(k,ie)) .gt. 0 ) then
	      i = ibhnd(k,ie)
              ieltv(i,ie) = -1
	    end if
          end do
        end do

!-------------------------------------------------------------------
! check structure
!-------------------------------------------------------------------

	call checkielt(nel,ieltv)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	end

!****************************************************************

        subroutine mkdxy(nkn,dxv,dyv)

! initializes dxv,dyv

        implicit none

! arguments
        integer nkn
	double precision dxv(nkn),dyv(nkn)
! local
	integer k,ie,ii,i
	integer ip,ip0,ip1

!-------------------------------------------------------------------
! loop over nodes
!------------------------------------------------------------------

        do k=1,nkn
	  dxv(k) = 0.
	  dyv(k) = 0.
        end do

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	end

        end module links
