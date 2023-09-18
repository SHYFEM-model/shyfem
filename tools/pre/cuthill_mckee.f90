!
! $Id: subcmv.f,v 1.8 2001/01/15 14:03:36 georg Exp $
!
! cuthill-mckee algorithm for optimization of bandwidth
!
! revision log :
!
! 20.03.1998	ggu	reorganized, automatic procedure introduced
! 05.06.1998	ggu	avoid write to terminal
!
!*******************************************************************
!--------------------------------------------------------------------------
        module cuthill_mckee
!--------------------------------------------------------------------------

        use utility
        use fem_util

!--------------------------------------------------------------------------
        contains
!--------------------------------------------------------------------------
	subroutine cmv(nkn,ngrddi,ipv,iphv,kphv,ng,iknot)

! cuthill-mckee algorithmus
!
! k..		neue knotennummern
! i..		alte knotennummern
! ng		grad des knoten
! iknot		kantenverzeichnis des knotens
!		...iknot(n,k) --> iknot(ngr*(k-1)+n)
!		...( n <= ng(k) )
! iphv		pointer fuer knotennummern : neu --> alt
! kphv		pointer fuer knotennummern : alt --> neu
! ikanf		anfangsknoten (erste stufe)
! ikmer		anfangsknoten fuer minimale bandbreite
! mmin		bandbreite fuer anfangsknoten ikmer
! knum		anzahl der bereits nummerierten knoten
! kanf,kend	neue knotennummern der alten stufe

        implicit none

	integer nkn,ngrddi
	integer ipv(nkn)
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)

        integer iwei
        integer jgrmin,jgrmax
        integer knum,m

! get switch iwei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iwei = 1

	do while( iwei .eq. 1 .or. iwei .eq. 2 )

	  call cmvwhat(iwei)

          if(iwei.eq.1) then				!grades
	    call getgrds(ngrddi,nkn,ng,jgrmin,jgrmax)
	    call cmgrade(nkn,ngrddi,ipv,iphv,kphv,ng,iknot,jgrmin,jgrmax)
	  else if(iwei.eq.2) then			!give first level
	    call getfstl(nkn,iphv,kphv,knum)
	    if( knum .gt. 0 ) then
              call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
	      write(6,*) 'node =',ipv(iphv(1)),'   mbw =',m
	    end if
	  else if(iwei.eq.3) then			!return old numbering
	    call ininum(nkn,iphv,kphv)
	  end if

	end do	!do while

	end

!********************************************************************

        subroutine cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)

! cm-algorithmus
!
! knum		anzahl der vorgegebenen knoten der ersten stufe
! kanf,kend	neue knotennummern der alten stufe
!		(...knum=kend-kanf+1)
! 		iphv,kphv muss fuer erste stufe schon
!		...besetzt sein (der rest null)
! m		enthaelt am ende gefundene bandbreite

        implicit none

        integer nkn,ngrddi,knum,m
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)

        integer kanf,kend,knold
        integer ngr,ka,inh,in
        integer nggmin,ngg,n,ia,mh

        ngr=ngrddi

	m=0
	kanf=1			!grenzen fuer
	kend=knum		!...erste stufe
        knold=0

        do while(knum.gt.knold)
          knold=knum
          do ka=kanf,kend
            ia=iphv(ka)   !knoten der alten stufe
            inh=1
            do while(inh.ne.0)  !inh=0 ==> zu ia gibt es
              inh=0     !...keinen unnummerierten
              nggmin=ngr+1    !...knoten mehr --> neues ia
              do n=1,ng(ia)   !suche knoten in der neuer
                in=iknot((ia-1)*ngr+n)  !...stufe der mit ia
                if(kphv(in).eq.0) then  !...verbunden und noch
                  ngg=ng(in)  !...nicht nummeriert
                  if(ngg.lt.nggmin) then  !...ist und
                    nggmin=ngg  !...kleinsten
                    inh=in  !...grad besitzt und
                  end if    !...und schreibe
                end if      !...ihn auf inh
              end do
              if(inh.ne.0) then !naechster zu
                knum=knum+1 !...nummerierender
                iphv(knum)=inh  !...knoten ist
                kphv(inh)=knum  !...inh und
                mh=iabs(ka-knum)!...die bandweite
                if(mh.gt.m) m=mh!...ist mh
              end if
            end do  !do while(inh.ne.0)
          end do  !ka
          kanf=kend+1   !neue stufe wird
          kend=knum   !...alte stufe
	end do	!do while

        return
	end

!****************************************************************

	subroutine cmvwhat(iwei)

! asks about what to do (CMV algorithm)

	implicit none

	integer iwei

	integer ianz
	double precision f(10)

	write(6,*) '0   save numbering and exit'
	write(6,*) '1   numbering for grades'
	write(6,*) '2   numbering starting from first level'
	write(6,*) '3   quit with numbering prior to routine call'

	iwei = -1
	do while( iwei .lt. 0 .or. iwei .gt. 3 )
	  iwei = 0
	  ianz = inquire_numbers('give number :',f)
	  if( ianz .gt. 0 ) iwei = iround(f(1))
	end do

	end

!****************************************************************

	subroutine getgrds(ngr,nkn,ng,jgrmin,jgrmax)

! gets min/max grades to use for a starting point minimization

	implicit none

	integer ngr,nkn
	integer ng(nkn)
	integer jgrmin,jgrmax

	double precision f(10)
	integer ianz,jz,n,j

	do n=1,ngr
            jz=0
            do j=1,nkn
              if(ng(j).eq.n) jz=jz+1
            end do

            if(jz.ne.0) then
              write(6,*) 'grade =',n,'   number of nodes =',jz
            end if
	end do

	ianz = -1
	do while( ianz .lt. 0 .or. ianz .gt. 2 )
	  ianz=inquire_numbers('give min/max grade (default 1/4) :',f)
	  if(ianz.eq.2) then
		jgrmin=f(1)
		jgrmax=f(2)
	  else if(ianz.eq.1) then
		jgrmin=f(1)
		jgrmax=f(1)
	  else if(ianz.eq.0) then
		jgrmin=1
		jgrmax=4
	  else
		write(6,*) 'erroneous input'
	  end if
	end do

	end

!****************************************************************

	subroutine getfstl(nkn,iphv,kphv,knum)

! gets first level interactively

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)
	integer knum

	integer ianz
	integer node,noden
	double precision f(10)

	knum=0
	call zernum(nkn,iphv,kphv)

	write(6,*)'give nodes of first level (<cr> to end)'

	ianz=inquire_numbers(' : ',f)
	do while(ianz.gt.0)
		if(ianz.gt.1) then
		  write(6,*)'only one node a time -> repeat'
		else
		  node = iround(f(1))
		  noden = ipint(node)
		  if(noden.ne.0) then
			knum=knum+1
			iphv(knum)=noden
			kphv(noden)=knum
		  else
			write(6,*)'invalid node number'
		  end if
		end if
		ianz=inquire_numbers(' : ',f)
	end do

	if(knum.eq.0) then	!no node
		write(6,*)'no valid node - old numbering reinstalled'
		call ininum(nkn,iphv,kphv)
	end if

	end

!****************************************************************

        subroutine nodnum(nkn,iphv,kphv,node,knum)

! initializes numbering of nodes for first node

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)
	integer node,knum

        integer i

        do i=1,nkn
          iphv(i) = 0
          kphv(i) = 0
        end do

        knum=1
        iphv(knum) = node
        kphv(node) = knum

	end

!****************************************************************

	subroutine cmgrade(nkn,ngrddi,ipv,iphv,kphv,ng,iknot,jgrmin,jgrmax)

! cuthill-mckee algorithm for grades

        implicit none

	integer nkn,ngrddi
	integer ipv(nkn)
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)
        integer jgrmin,jgrmax

        integer i
        integer mmin,ikmer,knum,m

	ikmer=0
	mmin=nkn
	m = 0

	write(6,*) 'Applying Cuthill-McKee algorithm...'
	!write(6,*) jgrmin,jgrmax
!        write(6,'(16x,a,7x,a,3x,a)') 'node','grade','bandwidth'

        do i=1,nkn
            if(ng(i).ge.jgrmin.and.ng(i).le.jgrmax) then
              call nodnum(nkn,iphv,kphv,i,knum)
              call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
	      !write(6,*) i,ipv(i),ng(i),m,mmin
	      call optest('inside cmgrade',nkn,ipv,iphv,kphv)
!              write(6,'(8x,3i12)') ipv(i),ng(i),m
              if(m.lt.mmin) then
                mmin=m
                ikmer=i
              end if
            end if
	end do

	if(ikmer.eq.0) then	!no valid node found
	    write(6,*)'no node for this grade -> '
	    write(6,*)'old numbering reinstalled'
	    call ininum(nkn,iphv,kphv)
        else			!repeat mimimal node
            call nodnum(nkn,iphv,kphv,ikmer,knum)
            call cmalg(nkn,ngrddi,knum,m,iphv,kphv,ng,iknot)
            write(6,*) 'node =',ipv(iphv(1)),'   mbw =',m
	end if

	end

!****************************************************************

!--------------------------------------------------------------------------
        end module cuthill_mckee
!--------------------------------------------------------------------------
