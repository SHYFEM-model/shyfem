c
c $Id: map_influence.f,v 1.2 2010-03-08 17:46:45 georg Exp $
c
c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 18.06.2009    ggu     re-written for map of influence
c 05.11.2010    ggu     write all conz to one file, bug fix in vert_aver()
c 16.12.2010    ggu     restructured to use volume file
c 25.02.2011    ggu     restructured for Isabella to sum over concentrations
c
c****************************************************************

	program nossum

c sums all concentrations in .com file to one concentration in .con file

	implicit none
	include 'param.h'
	
	integer nsdim
	parameter (nsdim=3)	!number of tracers

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
	real cv3(nlvdim,nkndim)
	real cv3_all(nlvdim,nkndim)
        
	integer nkn1,nkn2,nel1,nel2,nlv2
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	real hev(neldim)

	logical berror
	logical bvol
	integer nvol
        
c--------- local variables ----------------------

        real sum,rt,rnull
	real ptresh,ctresh
	real raux
	integer nin
	integer nunit,it,nvers,ivar,nvar,ierr
	integer iapini,ideffi
	integer nvarnew
	integer nout1,nout2,nout3,nout4,nout5
	integer i,nread,k,l,is,nwrite
	integer nstate
    
	integer ifem_choose_file,ifem_open_file

	character*6 namel(nsdim),nam
        
	integer icheck
	character*50 file

c---------------------------------------------------------------
c set important parameters
c---------------------------------------------------------------

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nvers=3

c---------------------------------------------------------------
c initializing variables
c---------------------------------------------------------------

	nread=0
	nwrite=0
	rnull=0.
       
	do l=1,nlvdim
          do k=1,nkn
	    cv3_all(l,k)=0.
          end do
        end do

c-------------------------------------------------------------
c open input file and read headers
c---------------------------------------------------------------

        nin = ifem_open_file('.nos','old')
        call rhnos(nin,nvers,nkndim,neldim,nlvdim,nkn1,nel1,nlv,nvar
     +                          ,ilhkv,hlv,hev,title)

	nstate = nvar
	!if( nstate .ne. nsdim ) goto 97

c---------------------------------------------------------------
c open output files
c---------------------------------------------------------------

	nvarnew = 1

        call open_nos_file('conz_all','new',nout3)
        call whnos(nout3,nvers,nkn,nel,nlv,nvarnew,ilhkv,hlv,hev,title)

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
	   cv3_all(l,k) = cv3_all(l,k) + cv3(l,k)
	 end do
	end do

	icheck = mod(nread,nstate)
	if( icheck .eq. 0 ) icheck = nstate
	if( icheck .ne. is ) goto 98

	if( is .eq. nstate ) then

c	   -------------------------------------------------------
c	   all state variables read for one time step -> elaborate
c	   -------------------------------------------------------

           call wrnos(nout3,it,30,nlvdim,ilhkv,cv3_all,ierr)

	   nwrite = nwrite + 1
	   do l=1,nlvdim
             do k=1,nkn
	       cv3_all(l,k)=0.
             end do
           end do

	end if

	goto 300

c---------------------------------------------------------------
c end of time loop
c---------------------------------------------------------------

  100   continue

        close(nin)
	close(nout3)

	write(6,*)
	write(6,*) nread,' records read'
	write(6,*) nwrite,' records written'
	write(6,*)
	write(6,*) 'data written to files conz_all.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   95   continue
        write(6,*) 'error parameters in fvl file : '
        write(6,*) 'nkn: ',nkn,nkn2
        write(6,*) 'nel: ',nel,nel2
        write(6,*) 'nlv: ',nlv,nlv2
        stop 'error stop nosaver: nkn,nel,nlv'
   96   continue
        write(6,*) 'error parameters in nos file: '
        write(6,*) 'nkn: ',nkn,nkn1
        write(6,*) 'nel: ',nel,nel1
        stop 'error stop nosaver: nkn,nel'
   97	continue
	write(6,*) 'nstate,nsdim: ',nstate,nsdim
	write(6,*) 'nstate and nsdim must be equal'
	stop 'error stop: nsdim'
   98	continue
	write(6,*) ivar,nread,nstate,icheck
	stop 'error stop: error in reading records...'
	end

c***************************************************************

	subroutine vert_aver(nstate,nlvdim,nkn,cv1,cv2d,vol3,ilhkv)

c vertical averaging

	implicit none

	integer nstate,nkn,nlvdim
	real cv1(nlvdim,nstate,1)
	real cv2d(nstate,1)
	real vol3(nlvdim,1)
	integer ilhkv(1)

	integer is,k,l,lmax
        double precision c,v
        double precision cctot,vvtot

	do is=1,nstate
	  do k=1,nkn
	    cctot = 0.
	    vvtot = 0.
	    lmax = ilhkv(k)
	    do l=1,lmax
	      c = cv1(l,is,k)
	      v = vol3(l,k)
	      cctot = cctot + c*v
	      vvtot = vvtot + v
	    end do
	    cctot = cctot / vvtot
	    cv2d(is,k) = cctot
	  end do
	end do

	end

c***************************************************************

	subroutine comp_map(nstate,nlvdim,nkn,pt,ct,cvv,valri)

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
