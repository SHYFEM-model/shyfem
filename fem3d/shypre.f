
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

c pre-processing routine
c
c revision log :
c
c 30.08.1988	ggu	(itief in common, no rdpara)
c 22.11.1988	ggu	(new version 3, itief passed as actual param)
c 24.11.1988	ggu	(sp13f., descrr...)
c 30.11.1988	ggu	(back to sp13u.)
c 31.07.1990	ggu	(open all files explicitly)
c 08.10.1994	ggu	(newly designed -> use subroutines)
c 09.10.1994	ggu	(read from .grd files)
c 16.03.1995	ggu	(double precision in clockw)
c 06.03.1996	ggu	renumber also iarv in renel
c 08.09.1997	ggu	introduce raux,neaux for compiler warnings
c 20.03.1998    ggu     automatic optimization of bandwidth introduced
c 08.05.1998    ggu     always process whole file (idepth = 0)
c 18.05.1998    ggu     always process depths elementwise
c 18.05.1998	ggu	dont ask for description anymore
c 17.10.2001    ggu     accept also grd files with some missing data
c 18.10.2005    ggu     some error messages slightly changed
c 06.04.2009    ggu     read param.h
c 24.04.2009    ggu     new call to rdgrd()
c 04.03.2011    ggu     new routine estimate_grade()
c 30.03.2011    ggu     new routine check_sidei(), text in optest()
c 15.07.2011    ggu     calls to ideffi substituted
c 15.11.2011    ggu     new routines for mixed depth (node and elem), hflag
c 09.03.2012    ggu     delete useless error messages, handle nkn/nel better
c 29.03.2012    ggu     use ngr1 to avoid too small dimension for ngr
c 04.10.2013    ggu     in optest better error handling
c 30.07.2015    ggu     vp renamed to shypre
c 17.11.2017    ggu     implement output switches (quiet,silent,etc..)
c 13.04.2018    ggu     accepts partition to write bas file with node partition
c
c notes :
c
c could eliminate scrambling of iknot ->
c       no knscr()
c       pass ngr1 to bandop
c       change cmv,rosen
c
c**********************************************************

c==========================================================
	module mod_shypre
c==========================================================

	logical, save :: bopti
	logical, save :: bauto
	logical, save :: bnoopti
	logical, save :: bmanual

	logical, save :: binfo
	logical, save :: bquiet
	logical, save :: bsilent

c==========================================================
	end module mod_shypre
c==========================================================

c**********************************************************

        program shypre

	use basin
	use mod_shypre

	implicit none

	character*80 name
	character*80 file
	character*80 errfil
	character*80 errtex
        character*80 descrg,descra,basnam,grdfile
        logical bstop,bwrite,bdebug,bww

        integer, save, allocatable :: ipdex(:), iedex(:)
        integer, save, allocatable :: iphv(:), kphv(:)
        integer, save, allocatable :: iphev(:), iaux(:), ierank(:)
	integer, save, allocatable :: neaux(:,:)
        real, save, allocatable :: raux(:)
        real, save, allocatable :: hev(:)
        real, save, allocatable :: hkv(:)
        integer, save, allocatable :: ng(:)
        integer, save, allocatable :: kvert(:,:)

	integer, save, allocatable :: iknot(:,:)

	integer ianz,idepth,itief
	integer ie,k
	integer nat,net,ner,nb1,nb2,nb3,nb4,nb9

	integer nk,ne,nl,nn
        integer nknh,nelh,nli,nco
	integer nknddi,nelddi
	integer ngr1
	integer nrec
	integer nlidim,nlndim
	integer nne,nnl
	real hflag

	integer idefbas,ichanm,ifileo

	data bstop /.false./
	data errfil /'errout.dat'/

	bdebug = .false.

        dcorbas=0.
        dirnbas=0.
        descrg=' '
        descrr=' '
        descra=' '
c
	nb1=1
	nb2=2
	nb3=3
	nb4=4
	nb9=9
	net=5
	nat=6
	ner=99		!errout

	hflag = -999.

	itief=0		!0=read by element  1=read by node

c--------------------------------------------------------
c get name of basin
c--------------------------------------------------------

	call shypre_init(grdfile)

c--------------------------------------------------------
c open error file
c--------------------------------------------------------

	ner=ifileo(ner,errfil,'form','new')
	if(ner.le.0) then
		write(6,*) 'Cannot open error file'
		stop 'error stop : shypre'
	end if

c--------------------------------------------------------
c set parameters
c--------------------------------------------------------

	bwrite = .not. bquiet
	bww = .not. bsilent

	call pardef
	basnam = grdfile
	call delete_extension(basnam,'.grd')
	call putfnm('basnam',grdfile)

c--------------------------------------------------------
c always process whole file
c--------------------------------------------------------

	idepth = 0

c--------------------------------------------------------
c read grid
c--------------------------------------------------------

	if( bww ) write(6,*) 'grdfile: ',trim(grdfile)
	call grd_read(grdfile)

	call grd_get_params(nk,ne,nl,nne,nnl)
	if( bwrite ) write(6,*) 'grid info: ',nk,ne,nl

	if( nk == 0 .or. ne == 0 ) then
	  write(6,*) 'nk,ne: ',nk,ne
	  stop 'error stop shypre: no nodes or elements in basin'
	end if

	call grd_to_basin

c--------------------------------------------------------
c allocate arrays
c--------------------------------------------------------

	allocate(ipdex(nk), iedex(ne))
	allocate(iphv(nk), kphv(nk))
	allocate(iphev(ne), iaux(ne))
	allocate(neaux(3,ne))
	allocate(raux(ne))
	allocate(hev(ne))
	allocate(hkv(nk))
	allocate(ierank(ne))
	allocate(ng(nk))
	allocate(kvert(2,nk))

c--------------------------------------------------------
c handle depth
c--------------------------------------------------------

	call grd_get_depth(nk,ne,hkv,hev)

	nknh = 0
	do k=1,nkn
	  if( hkv(k) .ne. hflag ) nknh = nknh + 1
	end do

	nelh = 0
	do k=1,nel
	  if( hev(k) .ne. hflag ) nelh = nelh + 1
	end do

	if( bww ) then
          write(6,*) 'nkn,nel   : ',nkn,nel
          write(6,*) 'nknh,nelh : ',nknh,nelh
          !write(6,*) 'nli,nco   : ',nli,nco
	end if

        if(nkn.le.0 .or. nel.le.0) then
          write(ner,*) ' Nothing to process'
	  goto 99999
        end if

	itief=0
	if( nel.eq.nelh .and. nkn.eq.nknh ) then
	 if( bwrite ) then
	  write(nat,*) ' Can process depth node or elementwise.'
          write(nat,*) ' ...depths are processed elementwise'
	 end if
	else if(nel.eq.nelh) then
	 if( bwrite ) then
          write(nat,*) ' ...depths are processed elementwise'
	 end if
	else if(nkn.eq.nknh) then
	  itief=1
	 if( bwrite ) then
          write(nat,*) ' ...depths are processed nodewise'
	 end if
	else if(nknh.eq.0.and.nelh.eq.0) then
	 if( bwrite ) then
	  write(nat,*) ' No depth data read. Process anyway'
	 end if
	else
	  itief=2
	 if( bwrite ) then
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	  write(nat,*) ' Mixed data source for depth. Process anyway'
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	 end if
	end if

	if( binfo ) stop

c--------------------------------------------------------
c open files
c--------------------------------------------------------

	nb2=idefbas(basnam,'new')
        if(nb2.le.0) stop

c--------------------------------------------------------
c start processing
c--------------------------------------------------------

	nknddi = nkn
	nelddi = nel

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest('end read',nelddi,nkn,nel,nen3v)
	end if

        bstop=.false.

	if( bwrite ) then
        write(nat,*) ' ...making external index'
	end if

        call isort(nkn,ipv,ipdex)
        call isort(nel,ipev,iedex)

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of node numbers'
	end if

        call uniqn(nkn,ipv,ipdex,ner,bstop)
	if(bstop) goto 99909

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of element numbers'
	end if

        call uniqn(nel,ipev,iedex,ner,bstop)

	if( bdebug ) then
	call gtest('end uniqn',nelddi,nkn,nel,nen3v)
	end if

	if( bwrite ) then
        write(nat,*) ' ...controlling uniqueness of elements'
	end if

        call uniqe(nel,nen3v,iaux,iphev,ipev,ner,bstop)
        if(bstop) goto 99915

	if( bwrite ) then
	write(nat,*) ' ...changing extern with intern node numbers'
	end if

        !call chexin(nkn,nel,nen3v,ipv,ipdex,ner,bstop)
	!if(bstop) goto 99920

	if( bwrite ) then
        write(nat,*) ' ...controlling node numbers'
	end if

        call needn(nkn,nel,nen3v,ipv,iaux,ner,bstop)
	if(bstop) goto 99918

	if( bwrite ) then
	write(nat,*) ' ...testing sense of nodes in index'
	end if

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest('end sense',nelddi,nkn,nel,nen3v)
	end if

        call clockw(nkn,nel,nen3v,ipev,xgv,ygv,ner,bstop)
	if(bstop) goto 99930

	if( bwrite ) then
	write(nat,*) ' ...setting up side index'
	end if

        call estimate_grade(nkn,nel,nen3v,ng,ngr1)
	allocate(iknot(ngr1,nk))

        call sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr)
        call check_sidei(nkn,nel,nen3v,ipv,ng,iknot,ngr1,iaux)
	call knscr(nkn,ngr,ngr1,iknot)

	if( bww ) write(nat,*) 'Maximum grade of nodes is ',ngr

c--------------------------------------------------------
c bandwidth optimization
c--------------------------------------------------------

	if( bdebug ) then
	call gtest('bandwidth',nelddi,nkn,nel,nen3v)
	end if

        if( bwrite ) write(nat,*) ' ...optimizing band width'

        call bandw(nel,nen3v,mbw)

	if( bdebug ) then
	call gtest('bandwidth 1',nelddi,nkn,nel,nen3v)
	end if
	if( bww ) write(nat,*) 'Bandwidth is ',mbw

        call bandop(nkn,ngr,ipv,iphv,kphv,ng,iknot,kvert
     +			,bopti,bauto,bww)

	if( bdebug ) then
	call gtest('bandwidth 2',nelddi,nkn,nel,nen3v)
	end if

        if(bopti) then
          call bandex(nkn,nel,nen3v,kphv,ipv,iarnv,xgv,ygv,hkv) !kphv is n-rank
	  if( bdebug ) then
	    call gtest('bandwidth 3',nelddi,nkn,nel,nen3v)
	  end if
          call bandw(nel,nen3v,mbw)
	  if( bdebug ) then
	    call gtest('bandwidth 4',nelddi,nkn,nel,nen3v)
	  end if
	  if( bww ) write(nat,*) 'Optimized bandwidth is ',mbw
	end if

c--------------------------------------------------------
c renumber elements
c--------------------------------------------------------

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest('end bandwidth',nelddi,nkn,nel,nen3v)
	end if

	if( bwrite ) write(nat,*) ' ...renumbering elements'

        call renel(nel,nen3v,ipev,iarv,hev,ierank)	!ierank is element rank

c--------------------------------------------------------
c save pointers for depth
c--------------------------------------------------------

        if( bwrite ) write(nat,*) ' ...saving pointers'

	do k=1,nkn
	  iphv(k)=ipv(k)
	end do

	do ie=1,nel
	  iphev(ie)=ipev(ie)
	end do

	call ketest(nel,nen3v)
	if( bdebug ) then
	call gtest('write',nelddi,nkn,nel,nen3v)
	end if

c--------------------------------------------------------
c process depths
c--------------------------------------------------------

        if( bwrite ) write(nat,*) ' ...processing depths'

	call init_hm3v(nel,hm3v,hflag)

	if(itief.eq.0) then
          if( bwrite ) write(nat,*) ' ...................elementwise'
          call helem(nel,nelh,iphev,iaux,iedex,ipev
     +				,hm3v,hev,ner,bstop)
	else if( itief .eq. 1 ) then
          if( bwrite ) write(nat,*) ' ......................nodewise'
          call hnode(nkn,nel,nknh,nen3v,iphv,iaux,ipdex,ipv
     +                      ,hm3v,hkv,ner,bstop)
	else
          if( bwrite ) write(nat,*) ' .............first elementwise'
          call helem(nel,nelh,iphev,iaux,iedex,ipev
     +				,hm3v,hev,ner,bstop)
          if( bwrite ) write(nat,*) ' .................then nodewise'
          call hnode(nkn,nel,nknh,nen3v,iphv,iaux,ipdex,ipv
     +                      ,hm3v,hkv,ner,bstop)
	end if

	call check_hm3v(nel,hm3v,hflag)

	if( bstop ) then
	  write(6,*) '*** error in processing depth'
	end if

	call ketest(nel,nen3v)

	descrr=descrg

c--------------------------------------------------------
c process partitions
c--------------------------------------------------------

	call handle_partition(nkn,nel,kphv,ierank)

c--------------------------------------------------------
c write to file
c--------------------------------------------------------

	if( bwrite ) write(nat,*) ' ...writing file ',nb2

	call basin_write(nb2)

	close(nb2)

	if( bww ) call bas_info

c--------------------------------------------------------
c end of routine
c--------------------------------------------------------

	stop
99900	write(nat,*)' (00) error in dimension declaration'
	goto 99
99909	write(nat,*)' (09) error : no unique definition of nodes'
	goto 99
99915	write(nat,*)' (15) error : no unique definition of elements'
	goto 99
99918	write(nat,*)' (18) error reading file 1'
	goto 99
99920	write(nat,*)' (20) error reading file 1'
	goto 99
99930	write(nat,*)' (30) error in element index'
	goto 99
99990	continue
	write(6,*) 'nknddi = ',nknddi,'    nkn = ',nkn
	write(6,*) 'nelddi = ',nelddi,'    nel = ',nel
	write(nat,*)' (90) error reading grid file'
	write(nat,*)' plaese have a look at the error messages'
	write(nat,*)' if there are errors in the file please adjust'
	write(nat,*)' if there are dimension errors please adjust'
	write(nat,*)'   dimensions, re-compile and re-run the command'
	goto 99
99999	write(nat,*)' (99) generic error'
	goto 99
c
   99	continue
c
	nrec = 0
	rewind(ner)
   98	read(ner,6000,err=97,end=97) errtex
	write(nat,*) errtex
	nrec = nrec + 1
	goto 98
   97	continue
c
	if( nrec .gt. 0 ) then
	  write(nat,*)
	  write(nat,*) ' Error messages have been written to file :'
	  write(nat,*) errfil
	  write(nat,*)
	end if

	stop ' error stop : shypre'
 6000	format(a)
 	end

c*****************************************************************

        subroutine uniqn(n,iv,index,ner,bstop)

c controlls uniqueness of numbers in iv

        implicit none

        integer n,ner
        logical bstop
        integer iv(n)
        integer index(n)

        integer k1,k2,i

        k1=iv(index(1))
        do i=2,n
          k2=iv(index(i))
          if(k1.eq.k2) then
            write(ner,*)' warning : ',k1,' not unique'
            bstop=.true.
          end if
          k1=k2
	end do

        return
        end

c*****************************************************************

        subroutine uniqe(nel,nen3v,iaux,index,ipev,ner,bstop)

c controlls uniqueness of elements (shell for equale)

        implicit none

        integer nel,ner
        logical bstop
        integer nen3v(3,nel)
        integer iaux(nel),index(nel+1)  !bug lahey
        integer ipev(nel)

        integer ie,ii,ka,km,k1,k
        integer ie1,ie2
        integer isum

        isum=0

        do ie=1,nel
          ka=0
          km=1
          k1=1
          do ii=1,3
            k=nen3v(ii,ie)
            ka=ka+k
            km=km*(mod(k,97)+1)
            k1=k1*(mod(k+1,97)+1)
          end do
          iaux(ie)=ka+km+k1
        end do

        call isort(nel,iaux,index)

        ie1=1
        do while(ie1.le.nel)
          ie2=ie1+1
c          write(6,*) nel,ie1,ie2,index(ie1),index(ie2)
          !do while(ie2.le.nel.and.iaux(index(ie1)).eq.iaux(index(ie2)))
          do
	    if( ie2 > nel ) exit
	    if( iaux(index(ie1)) /= iaux(index(ie2)) ) exit
            ie2=ie2+1
          end do
          if(ie2-ie1.gt.1) then
            isum=isum+ie2-ie1-1
            call equale(nel,ie1,ie2-1,nen3v,index,ipev,ner,bstop)
          end if
          ie1=ie2
        end do

        !write(6,*) 'uniqe (isum) : ',isum

        return
        end

c*****************************************************************

        subroutine equale(nel,ip1,ip2,nen3v,index,ipev,ner,bstop)

c controlls uniqueness of listed elements

        implicit none

        integer nel
        integer ner
        integer ip1,ip2
        logical bstop
        integer nen3v(3,nel)
        integer index(nel)
        integer ipev(nel)

        integer ip,i,ie1,ie2
        integer i1,i2,i3,kn1,kn2,kn3

        do ip=ip1,ip2
          ie1=index(ip)
          kn1=nen3v(1,ie1)
          do i=ip+1,ip2
            ie2=index(i)
            do i1=1,3
              if(nen3v(i1,ie2).eq.kn1) then
                i2=mod(i1,3)+1
                i3=mod(i2,3)+1
                kn2=nen3v(2,ie1)
                kn3=nen3v(3,ie1)
                if(nen3v(i2,ie2).eq.kn2 .and.
     +                  nen3v(i3,ie2).eq.kn3) then
                  write(ner,*)' warning : element',ipev(ie1)
     +                  ,'  and',ipev(ie2),'  are identical'
                  bstop=.true.
                end if
              end if
            end do
          end do
	end do

        return
        end

c*****************************************************************

        subroutine needn(nkn,nel,nen3v,ipv,iaux,ner,bstop)

c controlls if all nodes are needed (use nen3v as one-dim array)

        implicit none

        integer nkn,nel,ner
        logical bstop
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer iaux(nkn)

        integer k,ie,ii

        do k=1,nkn
          iaux(k)=0
        end do

        do ie=1,nel
          do ii=1,3
            iaux(nen3v(ii,ie))=1
          end do
        end do

        do k=1,nkn
          if(iaux(k).eq.0) then
            write (ner,*)' warning : node ',ipv(k),' not needed'
            bstop=.true.
          end if
        end do

        return
        end

c*****************************************************************

        subroutine chexin(nkn,nel,nen3v,ipv,index,ner,bstop)

c changing extern with intern node numbers in element index

        implicit none

        integer nkn,nel,ner
        logical bstop
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer index(nkn)

        integer ie,ii,i,kn
        integer locate

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            i=locate(nkn,ipv,index,kn)
            if(i.le.0) then
              write(ner,*)' warning : node',kn,' not found'
              bstop=.true.
            else
              nen3v(ii,ie)=i
            end if
          end do
        end do

        return
        end

c*****************************************************************

        subroutine clockw(nkn,nel,nen3v,ipev,xgv,ygv,ner,bstop)

c test for anti-clockwise sense of nodes

        implicit none

        integer nkn,nel,ner
        logical bstop
        integer nen3v(3,nel)
        integer ipev(nel)
        real xgv(nkn),ygv(nkn)

        integer ie,ii,iii,k1,k2
        double precision a,x1,x2,y1,y2	!new 16.3.95

	do ie=1,nel
          a=0.
          do ii=1,3
            iii=mod(ii,3)+1
            k1=nen3v(ii,ie)
            k2=nen3v(iii,ie)
	    x1=xgv(k1)
	    x2=xgv(k2)
	    y1=ygv(k1)
	    y2=ygv(k2)
            a=a+x1*y2-x2*y1
          end do
          if(a.le.0.) then
            write(ner,*)' warning : nodes in element ',ipev(ie)
     +                     ,' are in clockwise sense'
            bstop=.true.
          end if
	end do

        return
        end

c*****************************************************************

        subroutine estimate_grade(nkn,nel,nen3v,ng,ngr)

c estimates ngr

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer ng(nkn)

        integer i,ii,ie,k,ngr

	do i=1,nkn
          ng(i)=0
	end do

	do ie=1,nel
          do ii=1,3
	    k = nen3v(ii,ie)
	    ng(k) = ng(k) + 1
	  end do
	end do

	ngr = 0
	do i=1,nkn
          ngr = max(ngr,ng(i))
	end do
	ngr = ngr + 1		!account for boundary nodes

	end

c*****************************************************************

        subroutine check_sidei(nkn,nel,nen3v,ipv,ng,iknot,ngr1,iaux)

c set up side index and find grade

        implicit none

        integer nkn,nel,ngr
        integer ngr1
        integer nen3v(3,nel)
        integer ipv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
        integer iaux(nkn)

        integer ii,ie,k
	integer igr,iel
	logical bstop

	bstop = .false.
	do k=1,nkn
	  iaux(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iaux(k) = iaux(k) + 1
	  end do
	end do

	do k=1,nkn
	  igr = ng(k)
	  iel = iaux(k)
	  if( igr .gt. iel + 1 .or. igr .lt. iel ) then
	    write(6,*) 'node not regular: ',ipv(k),igr,iel
	    bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop check_sidei: node not regular'

	end

c*****************************************************************

        subroutine sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr)

c set up side index and find grade ngr

        implicit none

        integer nkn,nel,ngr
        integer ngr1
        integer nen3v(3,nel)
        integer ng(nkn)
        integer iknot(ngr1,nkn)

        integer i,ii,iii,ie,k1,k2

	do i=1,nkn
	  do ii=1,ngr1
	    iknot(ii,i)=0
	  end do
          ng(i)=0
	end do

	do ie=1,nel
          do ii=1,3
            iii=mod(ii,3)+1
            k1=nen3v(ii,ie)
            k2=nen3v(iii,ie)
            do i=1,ng(k1)
              if(iknot(i,k1).eq.k2) goto 1    !already there
            end do
            ng(k1)=ng(k1)+1        !side not yet in index
	    if(ng(k1).gt.ngr1) goto 99
            iknot(ng(k1),k1)=k2
            ng(k2)=ng(k2)+1        !has to be added to index
	    if(ng(k2).gt.ngr1) goto 99
            iknot(ng(k2),k2)=k1
    1       continue
          end do
	end do

	ngr=0
	do i=1,nkn
          if(ng(i).gt.ngr) ngr=ng(i)
	end do

        return
   99	continue
	write(6,*) 'dimension of ngr1 is too low: ',ngr1
	stop 'error stop sidei: ngr1'
        end

c**********************************************************

	subroutine knscr(nkn,ngr,ngr1,iknot)

c scrambles side index for cmv and rosen

c -> please eliminate need for this

	implicit none

	integer nkn,ngr,ngr1
        integer iknot(ngr1*nkn)

	integer k,j

        do k=1,nkn
          do j=1,ngr
            iknot((k-1)*ngr+j)=iknot((k-1)*ngr1+j)
          end do
        end do

	return
	end

c**********************************************************

        subroutine bandw(nel,nen3v,mbw)

c determine bandwidth mbw

        implicit none

        integer nel,mbw
        integer nen3v(3,nel)

        integer ie,ii,iii,k,kk
        integer mh,mm

	mh=0

	do ie=1,nel
          do ii=1,3
            k=nen3v(ii,ie)
            do iii=ii+1,3
              kk=nen3v(iii,ie)
              mm=iabs(kk-k)
              if(mm.gt.mh) mh=mm
            end do
          end do
	end do

	mbw=mh

	return
	end

c**********************************************************

        subroutine bandop(nkn,ngr1,ipv,iphv,kphv,ng,iknot,kvert
     +				,bopti,bauto,bwrite)

c optimize band width

        implicit none

        integer nkn,ngr1
        logical bopti,bauto
        integer ipv(nkn)
        integer iphv(nkn),kphv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
        integer kvert(2,nkn)
	logical bwrite

	integer iantw

	if( bauto ) then
	  call ininum(nkn,iphv,kphv)
	  call optest('before optimization: ',nkn,ipv,iphv,kphv)
	  call cmgrade(nkn,ngr1,ipv,iphv,kphv,ng,iknot,1,4,bwrite)
	  call optest('after Cuthill McKee: ',nkn,ipv,iphv,kphv)
	  call revnum(nkn,iphv,kphv,bwrite)
	  call optest('after reversing nodes: ',nkn,ipv,iphv,kphv)
          call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
	  call optest('after Rosen: ',nkn,ipv,iphv,kphv)
	  return
	end if

        do while( bopti )
	  call ininum(nkn,iphv,kphv)

c	  call anneal(nkn,ngr1,kphv,ng,iknot,iphv,kvert)

          if( iantw(' Cuthill McKee algorithm ?') .gt. 0 ) then
            call cmv(nkn,ngr1,ipv,iphv,kphv,ng,iknot,bwrite)
	  end if

          if( iantw(' Reverse numbering of nodes ?') .gt. 0 ) then
	    call revnum(nkn,iphv,kphv,bwrite)
	  end if

          if( iantw(' Rosen algorithm ?') .gt. 0 ) then
            call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
	  end if

	  bopti = iantw(' Repeat optimization of bandwidth ?') .gt. 0
        end do

	bopti = .true.

        end

c**********************************************************

	subroutine revnum(nkn,iphv,kphv,bwrite)

c reverses numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)
	logical bwrite

	integer i

	if( bwrite ) write(6,*) 'Applying reverse algorithm...'

	do i=1,nkn
          kphv(i) = nkn+1-kphv(i)
          iphv(kphv(i)) = i
	end do

	end

c**********************************************************

	subroutine ininum(nkn,iphv,kphv)

c initializes numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)

	integer i

	do i=1,nkn
	  iphv(i) = i
	  kphv(i) = i
	end do

	end

c**********************************************************

	subroutine zernum(nkn,iphv,kphv)

c zeros numbering of nodes

	implicit none

	integer nkn
	integer iphv(nkn), kphv(nkn)

	integer i

	do i=1,nkn
	  iphv(i) = 0
	  kphv(i) = 0
	end do

	end

c**********************************************************

        subroutine bandex(nkn,nel,nen3v,kphv
     +                      ,ipv,iarnv,xgv,ygv,hkv)

c exchange nodes after optimization

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer kphv(nkn)
        integer ipv(nkn)
        integer iarnv(nkn)
        real xgv(nkn),ygv(nkn)
        real hkv(nkn)

        integer ie,ii
	integer neaux(3,nel)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ie)=kphv(neaux(ii,ie))
          end do
        end do

c       copy arrays with kphv as rank table

        call icopy(nkn,ipv,kphv)
        call icopy(nkn,iarnv,kphv)
        call rcopy(nkn,xgv,kphv)
        call rcopy(nkn,ygv,kphv)
        call rcopy(nkn,hkv,kphv)

        return
        end

c**********************************************************

        subroutine icopy(n,iv,irank)

c copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        integer iv(n)
        integer irank(n)

        integer i
	integer iauxv(n)

        do i=1,n
          iauxv(i)=iv(i)
        end do

        do i=1,n
          iv(irank(i))=iauxv(i)
        end do

        return
        end

c**********************************************************

        subroutine rcopy(n,rv,irank)

c copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        real rv(n)
        integer irank(n)

        integer i
	real rauxv(n)

        do i=1,n
          rauxv(i)=rv(i)
        end do

        do i=1,n
          rv(irank(i))=rauxv(i)
        end do

        return
        end

c**********************************************************

        subroutine renel(nel,nen3v,ipev,iarv,hev,ierank)

c renumbering of elements

c we construct iedex newly and use iaux,iedex as aux arrays
c iedex is also used as a real aux array for rcopy	-> changed
c neaux is probably he3v (real) used as an aux array	-> changed

        implicit none

        integer nel
        integer nen3v(3,nel)
        integer ipev(nel)
	integer iarv(nel)
        real hev(nel)
	integer ierank(nel)

        integer ie,ii
        integer iedex(nel)
        integer neaux(3,nel)
        integer ival(nel)

        do ie=1,nel
          ival(ie)=min(nen3v(1,ie),nen3v(2,ie),nen3v(3,ie))
	end do

        call isort(nel,ival,iedex)    !iedex is the index table
        call rank(nel,iedex,ierank)   !ierank is the rank table

        call icopy(nel,ipev,ierank)
        call icopy(nel,iarv,ierank)
        call rcopy(nel,hev,ierank)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ierank(ie))=neaux(ii,ie)
          end do
        end do

        return
        end

c**********************************************************

        subroutine rank(n,index,irank)

c builds rank table from index table

        implicit none

        integer n
        integer index(n),irank(n)

        integer i

        do i=1,n
          irank(index(i))=i
        end do

        return
        end

c**********************************************************

	subroutine init_hm3v(nel,hm3v,hinit)

	implicit none

	integer nel
	real hm3v(3,nel)
	real hinit

	integer ie,ii

	do ie=1,nel
	  do ii=1,3
	    hm3v(ii,ie) = hinit
	  end do
	end do

	end

c**********************************************************

	subroutine check_hm3v(nel,hm3v,hflag)

	implicit none

	integer nel
	real hm3v(3,nel)
	real hflag

	logical bmiss
	integer ie,ii,iflag
	real h

	iflag = 0

	do ie=1,nel
	  bmiss = .false.
	  do ii=1,3
	    h = hm3v(ii,ie)
	    if( h .eq. hflag ) then
	      iflag = iflag + 1
	      bmiss = .true.
	    end if
	  end do
	  if( bmiss ) write(6,*) ie,(hm3v(ii,ie),ii=1,3)
	end do

	if( iflag .gt. 0 ) then
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) 'flag found in depth: ',iflag,iflag/3
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	  write(6,*) '*******************************************'
	end if

	end

c**********************************************************

        subroutine helem(nel,nhd,iphev,iaux,iedex
     +                    ,ipev,hm3v,hev,ner,bstop)

c depth by elements

        implicit none

        integer nel,nhd
        integer ner
        logical bstop
        integer iaux(nel)
        integer iedex(nel)
        integer ipev(nel),iphev(nel)
        real hm3v(3,nel),hev(nel)

        integer ie,i,iel,ii
        integer locate

        do ie=1,nel
          iaux(ie)=0
        end do

	!write(6,*) 'helem: ',nel,nhd

        call isort(nel,ipev,iedex)

c it is : ipev(ie) == iphev(i)

        !do i=1,nhd
        do i=1,nel
          iel=iphev(i)
          ie=locate(nel,ipev,iedex,iel)
          if(ie.le.0) then
            write(ner,*)' warning : element',iel,' not found'
            bstop=.true.
          else
            if(iaux(ie).ne.0) then
              write(ner,*)' for element ',ipev(ie)
     +                        ,' depth data not unique'
              bstop=.true.
            else
              iaux(ie)=i
            end if
          end if
        end do

        do  ie=1,nel
          if(iaux(ie).eq.0) then
            write(ner,*)' for element ',ipev(ie)
     +                        ,' no depth data found'
            bstop=.true.
	  else
            do ii=1,3
              hm3v(ii,ie)=hev(iaux(ie))
            end do
	  end if
	end do

        return
        end

c**********************************************************

        subroutine hnode(nkn,nel,nhd,nen3v,iphv,iaux,ipdex,ipv
     +                      ,hm3v,hkv,ner,bstop)

c depth by nodes

        implicit none

        integer nkn,nel,nhd
        integer ner
        logical bstop
        integer iaux(nkn)
        integer ipdex(nkn)
        integer ipv(nkn),iphv(nkn)
        integer nen3v(3,nel)
        real hm3v(3,nel),hkv(nkn)

        integer ie,ii,i,k,kn
        integer locate
	real h,hflag

	hflag = -999.

        do k=1,nkn
          iaux(k)=0
        end do

        call isort(nkn,ipv,ipdex)

        do i=1,nhd
          kn=iphv(i)
          k=locate(nkn,ipv,ipdex,kn)
          if(k.le.0) then
            write(ner,*)' warning : node',kn,' not found'
            bstop=.true.
          else
            if(iaux(k).ne.0) then
              write(ner,*)' for node ',ipv(k)
     +                        ,' depth data not unique'
              bstop=.true.
            else
              iaux(k)=i
            end if
          end if
        end do

        do  k=1,nkn
          if(iaux(k).eq.0) then
            write(ner,*)' for node ',ipv(k)
     +                        ,' no depth data found'
            bstop=.true.
	  end if
	end do

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            h = hkv(iaux(kn))
            if( iaux(kn) .gt. 0 .and. h .ne. hflag ) then
              hm3v(ii,ie) = h
            end if
          end do
	end do

        return
        end

c**********************************************************

	subroutine gtest(text,nelddi,nkn,nel,nen3v)

	implicit none

	character*(*) text
	integer nelddi,nkn,nel
	integer nen3v(3,nelddi)

	integer ie,ii,k
	integer iii
	logical berror
	
	berror = .false.

	write(6,*) 'testing ... ',text
	write(6,*) nelddi,nkn,nel

	do ie=1,nel
	  do ii=1,3
	    iii = mod(ii,3) + 1
	    if( nen3v(ii,ie) .eq. nen3v(iii,ie) ) then
	      if( .not. berror ) then
		write(6,*) ie,(nen3v(iii,ie),iii=1,3)
		berror = .true.
	      end if
	    end if
	  end do
	end do

	if( berror ) then
	  write(6,*) 'testing has found errors...'
	end if
	end

c**********************************************************

	subroutine optest(text,nkn,ipv,iphv,kphv)

	implicit none

	character*(*) text
	integer nkn
	integer ipv(nkn)
	integer iphv(nkn)
	integer kphv(nkn)

	integer i

	do i=1,nkn
	  if( iphv(i) .gt. nkn ) goto 99
	  if( kphv(i) .gt. nkn ) goto 99
	  if( iphv(kphv(i)) .ne. i ) goto 99
	  if( kphv(iphv(i)) .ne. i ) goto 99
	end do

	return
   99	continue
	write(6,*) '*** Error in optest: '
	write(6,*) text
	write(6,*) 'error in pointers...'
	write(6,*) 'problem is close to following node...'
	write(6,*) 'node (intern/extern): ',i,ipv(i)
	write(6,*) i,nkn
	write(6,*) iphv(i),kphv(i)
	write(6,*) iphv(kphv(i)),kphv(iphv(i))
	write(6,*) 'maybe the domain is not connected...'
	stop 'error stop optest: pointers'
	end

c**********************************************************

	subroutine ketest(nel,nen3v)

c checks uniquness of nodes in elements

	implicit none

	integer nel
	integer nen3v(3,nel)

	integer ie,ii
	integer kn1,kn2,kn3

	do ie=1,nel
	  kn1 = nen3v(1,ie)
	  kn2 = nen3v(2,ie)
	  kn3 = nen3v(3,ie)
	  if( kn1 .eq. kn2 .or. kn1 .eq. kn3 .or. kn2 .eq. kn3 ) then
	    write(6,*) ie,kn1,kn2,kn3
	    stop 'error stop ketest: not unique nodes in element'
	  end if
	end do

	!write(77,*)
	!do ie=1,10
	!  write(77,*) ie,(nen3v(ii,ie),ii=1,3)
	!end do

	return
	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine handle_partition(nn,ne,knrank,ierank)

	use clo
	use grd
	use basin

	implicit none

	integer nn,ne
	integer knrank(nn)
	integer ierank(ne)

	logical bnepart,bgrd
	integer nnpart,nepart
	integer area_node(nn)
	integer area_elem(ne)

	integer i
	character*80 grdfile

        call clo_get_option('nepart',bnepart)
        call clo_get_option('partition',grdfile)
	bgrd = ( grdfile /= ' ' )

	if( .not. bnepart .and. .not. bgrd ) return

	if( bnepart .and. bgrd ) then
	  write(6,*) 'only one of -partition and -nepart can be given'
	  stop 'error stop handle_partition: options'
	end if

	if( bgrd ) then
	  write(6,*) 'reading partitioning file ',trim(grdfile)
	  call grd_read(grdfile)
	  if( nk_grd /= nn ) goto 99
	  if( ne_grd /= ne ) goto 99
	  area_node = ianv
	  area_elem = iaev
	else
	  area_node = iarnv
	  area_elem = iarv
	end if

	call renumber_partition(nn,area_node,nnpart)
	call renumber_partition(ne,area_elem,nepart)

	if( bgrd ) then
          call icopy(nn,area_node,knrank)
          call icopy(ne,area_elem,ierank)
	end if

	write(6,*) 'partitioning set: ',nnpart,nepart

	call basin_set_partition(nn,ne,nnpart,nepart,area_node,area_elem)

	return
   99	continue
	write(6,*) nk_grd,nn
	write(6,*) ne_grd,ne
	stop 'error stop handle_partition: incompatibility'
	end

c**********************************************************

	subroutine renumber_partition(n,area,npart)

	implicit none

	integer n
	integer area(n)
	integer npart

	integer i,ia,imax
	integer nmin,nmax
	integer, allocatable :: table_in(:),table_out(:)

	nmin = minval(area)
	nmax = maxval(area)

	if( nmin == nmax ) then		!no partition
	  npart = 0
	  area = 0
	  return
	end if

	if( nmax < 0 ) then
	  write(6,*) nmin,nmax
	  stop 'error stop renumber_partition: nmin,nmax'
	end if

	allocate(table_in(0:nmax),table_out(0:nmax))

	table_in = 0
	table_out = 0

	do i=1,n
	  ia = area(i)
	  table_in(ia) = table_in(ia) + 1
	end do

	imax = -1
	do ia=0,nmax
	  if( table_in(ia) > 0 ) then
	    imax = imax + 1
	    table_out(ia) = imax
	  end if
	end do
	npart = imax

	!write(6,*) 'nmax: ',nmax,npart
	!write(6,*) table_in
	!write(6,*) table_out

	!if( imax /= nmax ) then
	!  write(6,*) nmin,nmax
	!  stop 'error stop renumber_partition: internal error (1)'
	!end if

	do i=1,n
	  ia = area(i)
	  area(i) = table_out(ia)
	end do

	end

c**********************************************************
c**********************************************************
c**********************************************************

	subroutine shypre_init(grdfile)

	use clo
	use mod_shypre

	implicit none

	character*(*) grdfile

	call clo_init('shypre','grd-file','3.0')

        call clo_add_info('pre-processes grd file and create bas file')

	call clo_add_sep('general options')
	call clo_add_option('info',.false.,'only give info on grd file')
	call clo_add_option('quiet',.false.,'be quiet in execution')
	call clo_add_option('silent',.false.,'do not write anything')

	call clo_add_sep('optimization options')
	call clo_add_option('noopti',.false.,'do not optimize bandwidth')
        call clo_add_option('manual',.false.,'manual optimization')

	call clo_add_sep('options for partitioning')
	call clo_add_option('partition file',' '
     +		,'use file containing partitioning')
	call clo_add_option('nepart',.false.
     +		,'use node and elem type in file for partitioning')

	call clo_parse_options

        call clo_get_option('info',binfo)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)

        call clo_get_option('noopti',bnoopti)
        call clo_get_option('manual',bmanual)

	call clo_check_files(1)
	call clo_get_file(1,grdfile)

	bopti = .not. bnoopti
	bauto = .not. bmanual
	if( bsilent ) bquiet = .true.

	if( .not. bquiet ) then
          call shyfem_copyright('shypre - pre-processing of GRD grid')
	end if

	end

c**********************************************************

