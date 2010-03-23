c
c $Id: subini.f,v 1.7 2002/10/11 12:54:45 georg Exp $
c
c initialization and boundary input (2D)
c
c contents :
c
c subroutine volinl(kn,lmax,vol,conz,temp,salt)	inputs water into node
c subroutine voldist(dz,conz,temp,salt)		inputs water distributed
c
c subroutine rdrst(brest)			initializes values from restart
c subroutine iniuvz(const)			initializes variables
c
c revision log :
c
c 21.08.1998    ggu     copied from subn11.f
c 21.08.1998    ggu     xv eliminated
c 03.10.1998    ggu     initialize uov/vov during restart
c
c***************************************************************

	subroutine volinl(kn,lmax,vol,conz,temp,salt)

c inputs water into node (2d version)

	implicit none

	integer kn		!node where to input volume
	integer lmax		!maximum level for input (not used in 2d)
	real vol		!volume to input
	real conz		!concentration
	real temp		!temperature
	real salt		!salinity

        real coev(3,1)
        common /coev/coev
        real toev(3,1)
        common /toev/toev
        real soev(3,1)
        common /soev/soev

	real dz
	real vol1,vol2
	real mass1,mass2

c first input water volume

	call watvol(kn,vol1)
	call volnod(kn,vol,dz)
	call watvol(kn,vol2)
c	write(6,*) 'volinl1: ',vol1,vol2,vol2-vol1,vol,dz

c in dz is water level difference -> adjust concentrations...

	call connod(kn,dz,conz,coev)
	call connod(kn,dz,temp,toev)
	call connod(kn,dz,salt,soev)

	end

c*****************************************************************

	subroutine voldist(dz,conz,temp,salt)

c inputs water distributed over basin (2d version)

	implicit none

	real dz			!distributed water level rise
	real conz		!concentration
	real temp		!temperature
	real salt		!salinity

        real coev(3,1)
        common /coev/coev
        real toev(3,1)
        common /toev/toev
        real soev(3,1)
        common /soev/soev

c first input water

	call surel(0,dz)

c adjust concentrations...

	call conele(0,dz,conz,coev)
	call conele(0,dz,temp,toev)
	call conele(0,dz,salt,soev)

	end

c*****************************************************************

	subroutine rdrst(brest)

c reads and initializes values from restart

	implicit none

	logical brest		!on return true if restarted

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	real xv(3,1)
	common /xv/xv
        real uov(1),vov(1),unv(1),vnv(1)
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1),zeov(3,1)			!$$zeov
        common /zenv/zenv, /zeov/zeov
        real zov(1), znv(1)
        common /zov/zov, /znv/znv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v

	character*80 name,cdum
	integer nb,ier,nk,ne,nvers,idum
	integer k,ie,ii
	integer ifileo,rfout,rdout
	real dum

	brest = .false.

        call getfnm('restrt',name)

	if(name.eq.' ') return

c open file

	nb=ifileo(12,name,'unform','old')
	if(nb.le.0) goto 96

	write(6,*) 'reading restart-file :'
	write(6,*) name

c read header

	ier=rfout(nb,nvers,idum,idum,idum,idum,dum,dum,cdum)
	if(ier.ne.0) goto 98

c look for time record

	ier=0
	nk=nkn
	ne=nel
	it=itanf-1
	if(itanf.eq.-1) then
		do while(ier.eq.0)
		   ier=rdout(nb,nvers,it,nk,ne,xv,zenv,unv,vnv)
		end do
		if(ier.ne.-1) goto 98
		itanf=it
	else
		do while(ier.eq.0.and.it.lt.itanf)
		   ier=rdout(nb,nvers,it,nk,ne,xv,zenv,unv,vnv)
		end do
		if(ier.ne.0) goto 98
		if(it.gt.itanf) goto 98
	end if

c set some other parameters

	do k=1,nkn
	  up0v(k) = xv(1,k)
	  vp0v(k) = xv(2,k)
	  znv(k) = xv(3,k)
	  zov(k) = xv(3,k)
	end do

	do ie=1,nel			!$$zeov
	  do ii=1,3
	    zeov(ii,ie)=zenv(ii,ie)
	  end do
	  uov(ie) = unv(ie)
	  vov(ie) = vnv(ie)
	end do

c restart ok

	write(6,*) 'program restarted from it =',it
	write(6,*)
	close (nb)

	brest = .true.

	return
   96	continue
	write(6,*) 'error in opening restart file :'
	write(6,*) name
	stop 'error stop : rdrst'
   98	continue
	write(6,*) 'restart from it = ',itanf,'   not possible'
	stop 'error stop : rdrst'
	end

c*******************************************************************

	subroutine iniuvz(const)

c initializes variables

	implicit none

	real const	!constant z value to impose

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it

	integer nen3v(3,1)
	common /nen3v/nen3v
        real uov(1),vov(1),unv(1),vnv(1)
	common /uov/uov, /vov/vov, /unv/unv, /vnv/vnv
        real zenv(3,1),zeov(3,1)			!$$zeov
	common /zenv/zenv, /zeov/zeov
        real zov(1), znv(1)
        common /zov/zov, /znv/znv
        real up0v(1), vp0v(1)
        common /up0v/up0v, /vp0v/vp0v

	character*80 name
	integer k,ie,i,nb
	integer ifileo

        call getfnm('bound',name)

	if(name.eq.' ') then
		do k=1,nkn
		   znv(k)=const
		end do
	else
		nb=ifileo(10,name,'unform','old')
		if(nb.le.0) goto 81
		read(nb,err=82,end=82) (znv(k),k=1,nkn)
		close(nb)
		write(6,*) 'Initial water levels read from file : '
		write(6,*) name
	end if

	do k=1,nkn
	   up0v(k) = 0.
	   vp0v(k) = 0.
	   zov(k) = znv(k)
	end do

        do ie=1,nel
          do i=1,3
	    k = nen3v(i,ie)
            zenv(i,ie)=znv(k)
	    zeov(i,ie)=znv(k)
          end do
        end do

c new part - set u/v initial values in new arrays

        do ie=1,nel
          uov(ie)=0.
          unv(ie)=0.
          vov(ie)=0.
          vnv(ie)=0.
        end do

	return
   81	continue
	write(6,*) 'Error opening initial water level file :'
	write(6,*) name
	write(6,*) 'on unit ',nb
	stop 'error stop : iniuvz'
   82	continue
	write(6,*) 'Error reading from initial water level file'
	write(6,*) name
	stop 'error stop : iniuvz'
	end

c*****************************************************************
