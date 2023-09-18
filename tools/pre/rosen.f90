!
! $Id: subros.f,v 1.5 1998/06/18 11:21:44 georg Exp $
!
! Rosen bandwidth optimization algorithm
!
! revision log :
!
! 05.06.1998    ggu     avoid write to terminal
!
!*********************************************************************
!---------------------------------------------------------------------------
        module rosen
!---------------------------------------------------------------------------
        contains
!---------------------------------------------------------------------------

	subroutine rosen_alg(nkn,ngrddi,iphv,kphv,ng,iknot,kvert)
!
! rosen algorithmus
!
! k..		neue knotennummern
! i..		alte knotennummern
! ng		!grad des knotens
! iknot		!kantenverzeichnis des knotens
! 		...iknot(n,k) --> iknot(ngr*(k-1)+n)
!		...( n <= ng(k) )
! iphv		pointer fuer knotennummern : neu --> alt
! kphv		pointer fuer knotennummern : alt --> neu
! ipaar		indexpaare die die maximale bandbreite bestimmen
! nppar   anzahl der indexpaare auf ipaar (max. ndim)
! npaara	gesammtanzahl der indexpaare die die
!		...maximale bandbreite bestimmen
! khil		indexpaare die bei vertauschung die bandbreite
!		...unveraendert lassen
! nhil    anzahl der indexpaare auf khil (max. ndim)
! kvert		indexpaare die bereits miteinander
!		...vertauscht worden sind
! nvert		anzahl der indexpaare auf kvert (max. nkn)
!		...wird auf 0 gesetzt wenn es gelingt zwei
!		...indizes zu vertauschen bei denen die
!		...bandbreite erniedrigt wird
! lvert		=1 ==> es ist gelungen, ein indexpaar zu
!		...vertauschen. dabei ist nicht notwendigerweise
!		...die bandweite erniedrigt worden. das program
!		...wird beendet wenn in do-90-loop kein index-
!		...paar mehr vertauscht werden kann.
!
        integer ndim
        parameter (ndim=10)

	integer nkn,ngrddi
	integer iphv(nkn),kphv(nkn)
	integer ng(nkn),iknot(ngrddi*nkn)
	integer kvert(2,nkn)

        integer ipaar(2,ndim),khil(2,ndim)

        integer nhil,nvert,ngr
        integer m,npaar,npaara,lvert
        integer i,ii,kn1,kn2
        integer mh,ms,msloop
        integer ipp,ipa1,ipa2,kk,kmax
        integer kanf,k,mhhmax,mhh,mmin
        integer kmerk,k1,k2,kmin

        logical true
        integer msmax
        save true,msmax
        data true  / .true. /
        data msmax /  100   /
!
	write(6,*) 'Applying Rosen algorithm...'

	msloop=0
	nhil=0
	nvert=0
        ms=nkn+1
	ngr=ngrddi
	kmerk = 0
!
        do while (true)
!
! bestimmung der indexpaare die maximale
! ...bandbreite m haben und ablegen auf %%%%%%%%%%%%%%%%%%
! ...ipaar (max. ndim)
!
	m=0
	npaar=0
	npaara=0
	lvert=0
	do i=1,nkn
	kn1=kphv(i)
	do ii=1,ng(i)
	kn2=kphv(iknot(ngr*(i-1)+ii))
        mh=abs(kn2-kn1)
	if(mh.ge.m) then
		if(mh.eq.m) then
            if(npaar.lt.ndim) then
				npaar=npaar+1
				ipaar(1,npaar)=i
				ipaar(2,npaar)=iknot(ngr*(i-1)+ii)
			end if
			npaara=npaara+1
		else
			m=mh
			npaar=1
			npaara=1
			ipaar(1,1)=i
			ipaar(2,1)=iknot(ngr*(i-1)+ii)
		end if
	end if
	end do
	end do

!        write(6,'(a,i4,a,i4,a,i4,a,i4)')
!     +       ' m = '            ,m
!     +      ,' pairs = '        ,npaara
!     +      ,' yet exchanged = ',nvert
!     +      ,' tries = '        ,msloop

! new on may 93 t stop loop

        if(m.lt.ms) then
          ms=m
          msloop=msmax
        else if(msloop.eq.0) then
          npaar=0   !do not enter loop -> stop
        else
          msloop=msloop-1
        end if
!
! schleife ueber paare die erniedrigt werden muessen %%%%%%%%%%%%
!
	do ipp=1,npaar
!
	ipa1=ipaar(1,ipp)
	ipa2=ipaar(2,ipp)
	kn1=kphv(ipa1)			!kleinerer index
	kn2=kphv(ipa2)			!groesserer index
        if(kn1.gt.kn2) then
          call iswap(kn1,kn2)
          call iswap(ipa1,ipa2)
        end if
	mh=kn2-kn1			!bandbreite
	if(mh.lt.m) goto 90		!naechstes paar
	if(mh.gt.m) stop'error stop : m to big'
!
! groesseren index kn2 mit einem kleineren vertauschen
!
	kmax=kn1			!maximale knoten-
	do ii=1,ng(ipa2)		!...nummer bestimmen
          kk=kphv(iknot(ngr*(ipa2-1)+ii)) !...die mit kn2
          if(kk.gt.kmax) kmax=kk    !...verbunden ist
	end do
!
	kanf=kmax-m
        if(kanf.lt.1) kanf=1
	mmin=m				!minimale bandbreite
	do k=kn2-1,kanf,-1		!k ist index der mit
	ia=iphv(k)			!...kn2 vertauscht wird
	mhhmax=max(kmax-k,k-kn1)	!...und bandbreite
	if(mhhmax.gt.mmin) goto 10	!...mhhmax ergibt
	do ii=1,ng(ia)		!bestimmung der bandbreite fuer
	mhh=iabs(kn2-kphv(iknot(ngr*(ia-1)+ii))) !...restl. knoten
	if(mhh.eq.0) mhh=kn2-k	!kn2 war im knotenverzeichnis
	if(mhh.gt.mmin) goto 10
	if(mhh.gt.mhhmax) mhhmax=mhh
	end do
	if(mhhmax.lt.mmin) then		!kleinere bandbreite
		mmin=mhhmax		!...gefunden, merke
		kmerk=k			!...k und m
	else if(mhhmax.eq.m) then	!vertauschung laesst
          if(nhil.lt.ndim) then !...bandbreite
			nhil=nhil+1	!...gleich,
			khil(1,nhil)=k	!...trage in
			khil(2,nhil)=kn2!...nhil ein
		end if
	end if
   10	continue
	end do
!
	if(mmin.lt.m) then		!kleinere bandbreite gefunden
		nvert=0			!...--> vertausche
		call vertau(nkn,kmerk,kn2,iphv,kphv,lvert,nhil)
		goto 90
	end if
!
! kleineren index kn1 mit einem groesseren vertauschen
!
	kmin=kn2			!minimale knoten-
	do ii=1,ng(ipa1)		!...nummer bestimmen
	kk=kphv(iknot(ngr*(ipa1-1)+ii))	!...die mit kn1
	if(kk.lt.kmin) kmin=kk		!...verbunden ist
	end do
!
	kend=kmin+m
        if(kend.gt.nkn) kend=nkn
	mmin=m
	do k=kn1+1,kend		!k ist index der mit
	ia=iphv(k)			!...kn1 vertauscht wird
	mhhmax=max(k-kmin,kn2-k)
	if(mhhmax.gt.mmin) goto 11
	do ii=1,ng(ia)
	mhh=iabs(kn1-kphv(iknot(ngr*(ia-1)+ii)))
	if(mhh.eq.0) mhh=k-kn1	!kn1 war im knotenverzeichnis
	if(mhh.gt.mmin) goto 11
	if(mhh.gt.mhhmax) mhhmax=mhh
	end do
	if(mhhmax.lt.mmin) then
		mmin=mhhmax
		kmerk=k
	else if(mhhmax.eq.m) then
          if(nhil.lt.ndim) then
			nhil=nhil+1
			khil(1,nhil)=kn1
			khil(2,nhil)=k
		end if
	end if
   11	continue
	end do
!
	if(mmin.lt.m) then
		nvert=0
		call vertau(nkn,kn1,kmerk,iphv,kphv,lvert,nhil)
		goto 90
	end if
!
! bandbreite kann nicht erniedrigt werden
!
	if(nvert.ge.nkn) stop 'error stop :dimension kvert'
!
	do i=1,nhil		!ueberpruefen ob es noch
	do ii=1,nvert		!...indizes gibt auf nhil
	if(khil(1,i).eq.kvert(1,ii)) then!...die noch nicht
		if(khil(2,i).eq.kvert(2,ii)) goto 21	!...
	end if			!...vertauscht worden sind
	end do
	k1=khil(1,i)
	k2=khil(2,i)
	nvert=nvert+1
	kvert(1,nvert)=k1
	kvert(2,nvert)=k2
	call vertau(nkn,k1,k2,iphv,kphv,lvert,nhil)
	goto 90
   21	continue
	end do
!
! keine vertauschung moeglich
!
   90	continue
	end do
!
	if(lvert.eq.0) then	!keine vertauschung mehr moegl.
		write(6,*) 'minimal bandwidth found =',m
		write(6,*) 'total number of pairs =',npaara
		return
	end if
!
	end do	!do while
!
	end

!*****************************************************************

	subroutine vertau(nkn,k1,k2,iphv,kphv,lvert,nhil)

! exchanges indices

	implicit none

	integer nkn
	integer k1,k2
	integer iphv(nkn),kphv(nkn)
	integer lvert,nhil

	integer i1,i2

	lvert=1
	nhil=0
	i1=iphv(k1)
	i2=iphv(k2)
	iphv(k1)=i2
	iphv(k2)=i1
	kphv(i1)=k2
	kphv(i2)=k1

	return
	end

!******************************************************************

        subroutine iswap(a,b)

! swaps two integers

        integer a,b,c

        c=a
        a=b
        b=c

        return
	end
!---------------------------------------------------------------------------
        end module rosen
!---------------------------------------------------------------------------
