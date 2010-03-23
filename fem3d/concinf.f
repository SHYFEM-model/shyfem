c
c $Id: concinf.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c 18.11.1998    ggu     check dimensions with dimnos
c 06.04.1999    ggu     some cosmetic changes
c
c**************************************************************

	program nosinf

c reads nos file

	implicit none

	include 'param.h'

	integer ntdim
	parameter(ntdim=10625)

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv
        
	real cv(nkndim)
	real cv3(nlvdim,nkndim)

	real rval(ntdim,nkndim)

	integer ilhkv(nkndim),type(nkndim)
	real hlv(nlvdim)
	real hev(neldim),volel(neldim)

	integer nread,nin,nout
	integer nvers
	integer nlv,nvar
	integer ierr
	integer it,ivar,id
	integer l,k
	character*80 title
	real rnull
	real cmin,cmax
	real val0,tau
	real taumin,taumax
	real high
        real isave   
       
      	integer iapini,ideffi
                     
c--------------------------------------------------------------
        isave=0
	nread=0
	rnull=0.

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

	nin=ideffi('datdir','runnam','.nos','unform','old')
	if(nin.le.0) goto 100

        nvers=3
	call rfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
        if(ierr.ne.0) goto 100

        write(6,*) 'nvers    : ',nvers
        write(6,*) 'nkn,nel  : ',nkn,nel
        write(6,*) 'nlv,nvar : ',nlv,nvar
        write(6,*) 'title    : ',title

        call dimnos(nin,nkndim,neldim,nlvdim)

	call rsnos(nin,ilhkv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100
c-------------------------------------------------------------------
c determino volume singoli elementi
c--------------------------------------------------------------------
        
        call volume(nel,volel)
        open(65,file='masscons.dat',status='unknown',form='formatted')
        open(66,file='mass%con.dat',status='unknown',form='formatted')
c----------------------------------------------------------------------
       
        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)
	
        do while(.true.)

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
c	write(6,*) 'time : ',it,ivar

	do l=1,nlv
	  do k=1,nkn
	    cv(k)=cv3(l,k)
	  end do
	  call mimar(cv,nkn,cmin,cmax,rnull)
c         write(6,*) 'level,cmin,cmax : ',l,cmin,cmax
	end do
c-------------------------------------------------------------
c calcolo volume totale della sostanza per ogni time step
c----------------------------------------------------------------
        call conzvol(cv,nel,it) 
c----------------------------------------------------------------
c         write(86,*) it,cv(1000),cv(2000),cv(3000),cv(4000)
c         write(86,*) it,cv(7298),cv(8282),cv(9674),cv(9999),cv(4128)

	if( nread .gt. ntdim ) stop 'error stop ntdim'

	do k=1,nkn
	  rval(nread,k) = cv(k)
	end do

	end do	!do while

  100	continue

 	write(6,*)
 	write(6,*) nread,' records read'
 	write(6,*)

	high = 1.e+30
	taumin = high
	taumax = -high
c---------------------------------------------------------------------
c determina elementi interni la laguna su cui calcolare time res
c--------------------------------------------------------------------
c        nout=76 
c        open(nout,file='ndtype.dat',status='unknown',form='formatted')
c        do k=1,nkn
c        read(nout,*)type(k)
c        end do
c        close(nout)
	do k=1,nkn
	  type(k) = 0
	end do
     	
        val0 = 100.
	do k=1,nkn
         if(type(k).eq.0.)then
          cv(k)=0
          goto 66 
         end if
          call reg(nread,rval(1,k),val0,tau)
	  if( tau .gt. high ) then
            write(6,*) 'tau... ',k,tau
	    tau = 0.
	  else if( tau .lt. 1.e-5 ) then
            write(6,*) 'tau... ',k,tau
	    tau = 1000.
	  else
	    tau = 1. / tau
	  end if
	  taumax = max(taumax,tau)
	  taumin = min(taumin,tau)
	  cv(k) = tau / 24.
66       continue         
	end do

 	write(6,*) 'alpha min/max: ',taumin,taumax
 	write(6,*) 'tau   min/max: ',1./taumin,1./taumax

	nin = 77
        nvers=3
	it = 0
	ivar = 10
	open(nin,file='timeres.con',status='unknown',form='unformatted')
	call wfnos(nin,nvers,nkn,nel,nlv,nvar,title,ierr)
	call wsnos(nin,ilhkv,hlv,hev,ierr)
	call wrnos(nin,it,ivar,1,ilhkv,cv,ierr)
	close(nin)
        close(65)
        close(66)
	end

c***************************************************************
c subroutine rfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine wfnos(iunit,nvers,nkn,nel,nlv,nvar,title,ierr)
c subroutine rsnos(iunit,ilhkv,hlv,hev,ierr)
c subroutine wsnos(iunit,ilhkv,hlv,hev,ierr)
c subroutine rdnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)
c subroutine wrnos(iunit,it,ivar,nlvdim,ilhkv,c,ierr)
c***************************************************************

        subroutine mimar(xx,n,xmin,xmax,rnull)

c computes min/max of vector
c
c xx            vector
c n             dimension of vector
c xmin,xmax     min/max value in vector
c rnull		invalid value

        implicit none

        integer n,i,nmin
        real xx(n)
        real xmin,xmax,x,rnull

	do i=1,n
	  if(xx(i).ne.rnull) goto 1
	end do
    1	continue

	if(i.le.n) then
	  xmax=xx(i)
	  xmin=xx(i)
	else
	  xmax=rnull
	  xmin=rnull
	end if

	nmin=i+1

        do i=nmin,n
          x=xx(i)
	  if(x.ne.rnull) then
            if(x.gt.xmax) xmax=x
            if(x.lt.xmin) xmin=x
	  end if
        end do

        end

c***************************************************************

c*******************************************************

        subroutine reg(n,val,val0,alpha)

        implicit none

        integer n
        real val(1)
        real val0
        real alpha

        integer i
        real x
        double precision xt,t2

        xt = 0.
        t2 = 0.

        do i=1,n
          x = min(val0,val(i))
          x = log(val0/x)
          xt = xt + x * i
          t2 = t2 + i * i
        end do

        alpha = xt / t2

        end

c*******************************************************
c****************************************************
	subroutine volume(nel,volel)

	include 'param.h'
        
        integer nel
        real hm3v(3,neldim)
        common /hm3v/hm3v
        integer ie,ii 
        real dhel,prof,dpel
        real volel(neldim)
        real area,areatr
         
        do ie=1,nel
         area = areatr(ie)
         dhel=0
         do ii=1,3
          dpel=hm3v(ii,ie)
          dhel=dhel+dpel
         end do
         prof=dhel/3
         volel(ie)=area*prof
c        write(6,*)area,'***',prof,'****',volel(ie),ie
        end do 
        end

c*********************************************************
c***********************************************************
 	subroutine conzvol(cv,nel,it)
        
        integer k,ie,nel,it
        
        integer nen3v(3,1)
        common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v
	real zenv(3,1)
	common /zenv/zenv
        real cvk,spcon 
        real cv(1) 
	double precision cvs,scon,ssave
	save ssave

	integer icall
	save icall
	data icall /0/

        scon=0
        do ie=1,nel
	 area = areatr(ie)/3.
         cvs=0.
         do ii=1,3
          k=nen3v(ii,ie)
	  htot = hm3v(ii,ie) + zenv(ii,ie)
          cvs=cvs+cv(k)*htot*area 
         end do
	 scon = scon + cvs
        end do
        if(icall.eq.0) then
	  ssave = scon
	end if
        spcon=scon/ssave
        write(66,*)it,scon,spcon

	icall = 1

        end  
c************************************************************************
c**********************************************************************
        function areatr(ie)

c determination of area of element

c ie            number of element (internal)
c areatr        element area (return value)

        real areatr
        integer ie

        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,1)
        common /nen3v/nen3v

        real aj
        integer ii,i1,i2,k1,k2

        aj=0.
        do ii=1,3
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
          k1=nen3v(i1,ie)
          k2=nen3v(i2,ie)
          aj=aj+xgv(k1)*ygv(k2)-xgv(k2)*ygv(k1)
        end do

        areatr = aj / 2.

        end
      


