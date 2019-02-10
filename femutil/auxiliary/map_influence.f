
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c revision log :
c
c 18.11.1998    ggu     check dimensions with dimnos
c 24.02.1999    ggu     use n2int for node number translation
c 03.12.2001    ggu     cleaned up, hakata bay
c 18.06.2009    ggu     re-written for map of influence
c 05.11.2010    ggu     write all conz to one file, bug fix in vert_aver()
c 16.12.2010    ggu     restructured to use volume file
c 10.11.2011    ggu     call to init_volume() changed for hybrid levels
c
c****************************************************************

	program distribu_deb

c creates map of influence

	use mod_depth
	use basin

	implicit none
	include 'param.h'
	
	integer nsdim
	parameter (nsdim=3)	!number of tracers

c--------------------------------------------------
	integer nlv


c--------------------------------------------------

	character*80 title
        real cvv(nlvdim,nsdim,nkndim)
	real cvv_sum(nlvdim,nsdim,nkndim)
	real cv2d(nsdim,nkndim)
	real cv3(nlvdim,nkndim)
	real cv3_all(nlvdim,nkndim)
	real valri(nlvdim,nkndim)
	real valri2d(nkndim)
	real vol3(nlvdim,nkndim)
        
	integer nkn1,nkn2,nel1,nel2,nlv2
	integer ilhkv(nkndim)
	real hlv(nlvdim)
	integer ilhkv2(nkndim)
	real hlv2(nlvdim)
	real hev2(neldim)
	real hl(nlvdim)


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
	integer i,nread,k,l,is
	integer nstate
    
	integer ifem_choose_file,ifem_open_file

	character*6 namel(nsdim),nam
        
	integer icheck
	character*50 file

c---------------------------------------------------------------
c set important parameters
c---------------------------------------------------------------

	ptresh = 30.	!threshold on percentage
	ctresh = 0.	!threshold on concentration - 0 for everywhere
	ctresh = 0.001	!threshold on concentration - 0 for everywhere
	ctresh = 0.	!threshold on concentration - 0 for everywhere
	ctresh = 1.e-7	!threshold on concentration - 0 for everywhere
	ctresh = 10.	!threshold on concentration - 0 for everywhere
	ctresh = 1.	!threshold on concentration - 0 for everywhere
	ctresh = 100.	!threshold on concentration - 0 for everywhere

c---------------------------------------------------------------
c open simulation and basin
c---------------------------------------------------------------

	if(iapini(3,0,0,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nvers=3

c       ----------------------------------------------------------
c       file containing volumes
c       ----------------------------------------------------------

        nvol = ifem_choose_file('.fvl','old')
        bvol = nvol .gt. 0

        if( bvol ) then
          write(6,*) 'volume file opened... using it'
          call rhnos(nvol,nvers,nkndim,neldim,nlvdim,nkn2,nel2,nlv2,nvar
     +                          ,ilhkv2,hlv2,hev2,title)
        else
          write(6,*) 'cannot open volume file... doing without'
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

        call init_sigma_info(nlv,hlv)

        call init_volume(nlvdim,nkn,nel,nlv,nen3v,ilhkv,hlv,hev,hl,vol3)

	nstate = nvar
	if( nstate .ne. nsdim ) goto 97

c-----------------------------------------------------------------
c check compatibility
c-----------------------------------------------------------------

        if( nkn .ne. nkn1 ) goto 96
        if( nel .ne. nel1 ) goto 96

        if( bvol ) then
          if( nkn .ne. nkn2 ) goto 95
          if( nel .ne. nel2 ) goto 95
          if( nlv .ne. nlv2 ) goto 95
          call check_equal_i('ilhkv',nkn,ilhkv,ilhkv2)
          call check_equal_r('hlv',nlv,hlv,hlv2)
          call check_equal_r('hev',nel,hev,hev2)
        end if

c---------------------------------------------------------------
c open output files
c---------------------------------------------------------------

	nvarnew = 1

        call open_nos_file('maps','new',nout1)
        call whnos(nout1,nvers,nkn,nel,nlv,nvarnew,ilhkv,hlv,hev,title)

        call open_nos_file('maps2d','new',nout2)
        call whnos(nout2,nvers,nkn,nel,1,nvarnew,ilhkv,hlv,hev,title)

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
           cvv_sum(l,is,k) = cvv_sum(l,is,k) + cv3(l,k)
	   cvv(l,is,k) = cv3(l,k)
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

           if( bvol ) call get_volume(nvol,it,nlvdim,ilhkv,vol3)

	   call comp_map(nstate,nlvdim,nkn,ptresh,ctresh,cvv,valri)
           call wrnos(nout1,it,367,nlvdim,ilhkv,valri,ierr)

	   call vert_aver(nstate,nlvdim,nkn,cvv,cv2d,vol3,ilhkv)
	   call comp_map(nstate,1,nkn,ptresh,ctresh,cv2d,valri2d)
           call wrnos(nout2,it,367,1,ilhkv,valri2d,ierr)

           call wrnos(nout3,it,30,nlvdim,ilhkv,cv3_all,ierr)

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
	close(nout1)
	close(nout2)
	close(nout3)

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

        call open_nos_file('map_final','new',nout4)
        call whnos(nout4,nvers,nkn,nel,nlv,nvarnew,ilhkv,hlv,hev,title)
	call comp_map(nstate,nlvdim,nkn,ptresh,ctresh,cvv,valri)
        call wrnos(nout4,it,367,nlvdim,ilhkv,valri,ierr)
	close(nout4)
        
        call open_nos_file('map_final2d','new',nout5)
        call whnos(nout5,nvers,nkn,nel,1,nvarnew,ilhkv,hlv,hev,title)
	call vert_aver(nstate,nlvdim,nkn,cvv,cv2d,vol3,ilhkv)
	call comp_map(nstate,1,nkn,ptresh,ctresh,cv2d,valri2d)
        call wrnos(nout5,it,367,1,ilhkv,valri2d,ierr)
	close(nout5)

c-------------------------------------------------------------
c final messages
c--------------------------------------------------------------

        if( .not. bvol ) then
          write(6,*)
          write(6,*) 'no volume file found: average done without'
        end if

	write(6,*)
	write(6,*) 'final data written to files map_final(2d).nos'
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
