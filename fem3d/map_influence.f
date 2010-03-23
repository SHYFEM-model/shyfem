c
c $Id: map_influence.f,v 1.2 2010-03-08 17:46:45 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 18.06.2009    ggu     re-written for map of influence
c 01.03.2010    ggu     adjusted and cleaned up
c
c****************************************************************

	program distribu_deb

c creates map of influence

	implicit none
	include 'param.h'
	
	integer nsdim
	parameter (nsdim=9)	!number of tracers

c--------------------------------------------------
        character*80 descrr
        common /descrr/descrr
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw,nlv
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        integer nen3v(3,neldim)
        integer ipv(nkndim), ipev(neldim)
        integer iarv(neldim)

        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v
        common /nen3v/nen3v
        common /ipv/ipv, /ipev/ipev
        common /iarv/iarv
c--------------------------------------------------

	character*80 title
        real cvv(nlvdim,nsdim,nkndim)
	real cvv_sum(nlvdim,nsdim,nkndim)
	real cv2d(nsdim,nkndim)
	real cv3(nlvdim,nkndim)
	real valri(nlvdim,nkndim)
	real valri2d(nkndim)
        
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
        
c--------- local variables ----------------------

        real sum,rt,rnull
	real ptresh,ctresh
	real raux
	integer nin
	integer nunit,nout,it,nvers,ivar,nvar,ierr
	integer iapini,ideffi
	integer nout2,nvarnew
	integer i,nread,k,l,is
	integer nstate
    
	character*6 namel(nsdim),nam
        
	integer icheck
	character*50 file

c---------------------------------------------------------------
c set important parameters
c---------------------------------------------------------------

	ptresh = 30.	!threshold on percentage
	ctresh = 100.	!threshold on concentration - 0 for everywhere
	ctresh = 0.	!threshold on concentration - 0 for everywhere
	ctresh = 0.001	!threshold on concentration - 0 for everywhere
	ctresh = 0.	!threshold on concentration - 0 for everywhere
	ctresh = 1.e-7	!threshold on concentration - 0 for everywhere

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c---------------------------------------------------------------
c initializing variables
c---------------------------------------------------------------

	nread=0
	rnull=0.
       
	do l=1,nlvdim
          do is=1,nsdim
            do k=1,nkn
	      cvv(l,is,k)=0.
	      cvv_sum(l,is,k)=0.
            end do
          end do
	end do

c-------------------------------------------------------------
c open input file and read headers
c---------------------------------------------------------------

	nin=ideffi('datdir','runnam','.nos','unform','old')

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

	write(6,*) 'Available levels: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

	nstate = nvar
	if( nstate .ne. nsdim ) goto 97

c---------------------------------------------------------------
c open output files
c---------------------------------------------------------------

        nvers = 3
	nvarnew = 1

	nout = 67
        open(nout,file='maps.nos',status='unknown'
     +				,form='unformatted')

	call wfnos(nout,nvers,nkn,nel,nlv,nvarnew,title,ierr)
        call wsnos(nout,ilhkv,hlv,hev,ierr)

	nout2 = 68
        open(nout2,file='maps2d.nos',status='unknown'
     +				,form='unformatted')

	call wfnos(nout2,nvers,nkn,nel,1,nvarnew,title,ierr)
        call wsnos(nout2,ilhkv,hlv,hev,ierr)
	
c	nout3 = 69
c       open(nout3,file='conz2d.nos',status='unknown'
c     +				,form='unformatted')

c	call wfnos(nout3,nvers,nkn,nel,1,nvar,title,ierr)
c       call wsnos(nout3,ilhkv,hlv,hev,ierr)
	
c---------------------------------------------------------------
c time loop
c---------------------------------------------------------------

  300   continue

	call rdnos(nin,it,ivar,nlvdim,ilhkv,cv3,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	write(6,*) 'time : ',it,ivar
	is = ivar - 30
	
	do l=1,nlv
	 do k=1,nkn
           cvv_sum(l,is,k) = cvv_sum(l,is,k) + cv3(l,k)
	   cvv(l,is,k) = cv3(l,k)
	 end do
	end do

	icheck = mod(nread,nstate)
	if( icheck .eq. 0 ) icheck = nstate
	if( icheck .ne. is ) ) goto 98

	if( is .eq. nstate ) then

c	   -------------------------------------------------------
c	   all state variables read for one time step -> elaborate
c	   -------------------------------------------------------

	   call comp_conz(nstate,nlvdim,nkn,ptresh,ctresh,cvv,valri)
           call wrnos(nout,it,367,nlvdim,ilhkv,valri,ierr)
	   call vert_aver(nstate,nlvdim,nkn,cvv,cv2d,ilhkv,hlv,hev)
	   call comp_conz(nstate,1,nkn,ptresh,ctresh,cv2d,valri2d)
           call wrnos(nout2,it,367,1,ilhkv,valri2d,ierr)

	end if

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100   continue

        close(nin)
	close(nout)
	close(nout2)

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*)
	write(6,*) 'data written to files maps(2d).nos'
	write(6,*)

c---------------------------------------------------------------
c compute average concentration over whole simulation
c---------------------------------------------------------------

	raux = 0.
	if( nread .gt. 0 ) raux = float(nstate) / float(nread)

	do is=1,nstate
	  do k=1,nkn
	    do l=1,nlvdim
              cvv(l,is,k) = cvv_sum(l,is,k) * raux
	    end do
	  end do
	end do

c---------------------------------------------------------------
c write final accumulated maps
c---------------------------------------------------------------

	nout=61
        nvers=3

        open(nout,file='map_final.nos',status='unknown'
     +				,form='unformatted')

	call wfnos(nout,nvers,nkn,nel,nlv,nvar,title,ierr)
        call wsnos(nout,ilhkv,hlv,hev,ierr)

	call comp_conz(nstate,nlvdim,nkn,ptresh,ctresh,cvv,valri)
        call wrnos(nout,it,367,nlvdim,ilhkv,valri,ierr)
	close(nout)
        
	call vert_aver(nstate,nlvdim,nkn,cvv,cv2d,ilhkv,hlv,hev)
	call comp_conz(nstate,1,nkn,ptresh,ctresh,cv2d,valri2d)

	nout2 = 62
        open(nout2,file='map_final2d.nos',status='unknown'
     +				,form='unformatted')

	call wfnos(nout2,nvers,nkn,nel,1,nvarnew,title,ierr)
        call wsnos(nout2,ilhkv,hlv,hev,ierr)
        call wrnos(nout2,it,367,1,ilhkv,valri2d,ierr)
	close(nout2)

c-------------------------------------------------------------
c final messages
c--------------------------------------------------------------

	write(6,*)
	write(6,*) 'final data written to files map_final(2d).nos'
	write(6,*)
	
c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   97	continue
	write(6,*) 'nstate,nsdim: ',nstate,nsdim
	write(6,*) 'nstate and nsdim must be equal'
	stop 'error stop: nsdim'
   98	continue
	write(6,*) ivar,nread,nstate,icheck
	stop 'error stop: error in reading records...'
	end

c***************************************************************

	subroutine vert_aver(nstate,nlvdim,nkn,cv1,cv2d,ilhkv,hlv,hev)

c vertical averaging

	implicit none

	integer nstate,nkn,nlvdim
	real cv1(nlvdim,nstate,1)
	real cv2d(nlvdim,1)
	integer ilhkv(1)
	real hlv(1),hev(1)

	integer is,k,l,lmax
	real conz
	real h,htop,hbot,htot

	do is=1,nstate
	  do k=1,nkn
	    lmax = ilhkv(k)
	    conz = 0.
	    htop = 0.
	    htot = 0.
	    do l=1,lmax
	      hbot = hlv(l)
	      h = hbot - htop
	      htot = htot + h
	      conz = conz + h * cv1(l,is,k)
	      htop = hbot
	    end do
	    cv2d(is,k) = 0.
	    if( htot .gt. 0. ) cv2d(is,k) = conz / htot
	  end do
	end do

	end

c***************************************************************

	subroutine comp_conz(nstate,nlvdim,nkn,pt,ct,cvv,valri)

c compute dominant discharge and put index in valri

	implicit none

	integer nstate,nkn,nlvdim
	real pt,ct
	real cvv(nlvdim,nstate,1)
	real valri(nlvdim,1)

	integer is,k,ismax,l
	real conz, pconz
	real sum,rmax
	real ctresh,ptresh

	ptresh = pt	!threshold on percentage
	ctresh = ct	!threshold on concentration - 0 for everywhere

	do l=1,nlvdim
	  do k=1,nkn
		sum = 0.
		rmax = 0.
		ismax = 0
        	do is=1,nstate
		   conz = cvv(l,is,k)
                   sum = sum + conz
		   if( conz .gt. rmax ) then
			ismax = is
			rmax = conz
		   end if
                end do

		conz = 0.
		if( ismax .gt. 0 ) conz = cvv(l,ismax,k)
		pconz = 0.
		if( sum .gt. 0. ) pconz = (conz/sum)*100

		valri(l,k) = 0.
		if( conz .gt. ctresh ) then
                  if( pconz .gt. ptresh ) then
		    valri(l,k) = ismax
 	          end if
		end if
	  end do
	end do

	end

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
