c
c $Id: vp.f,v 1.11 2009-05-21 09:24:00 georg Exp $
c
c pre-processing routine
c
c revision log :
c
c revised on  30.08.88  by ggu  (itief in common, no rdpara)
c revised on  22.11.88  by ggu  (new version 3, itief passed as actual param)
c revised on  24.11.88  by ggu  (sp13f., descrr...)
c revised on  30.11.88  by ggu  (back to sp13u.)
c revised on  31.07.90  by ggu  (open all files explicitly)
c revised on  08.10.94  by ggu  (newly designed -> use subroutines)
c revised on  09.10.94  by ggu  (read from .grd files)
c revised on  16.03.95  by ggu  (double precision in clockw)
c revised on  06.03.96  by ggu  renumber also iarv in renel
c revised on  08.09.97  by ggu  introduce raux,neaux for compiler warnings
c 20.03.1998    ggu     automatic optimization of bandwidth introduced
c 08.05.1998    ggu     always process whole file (idepth = 0)
c 18.05.1998    ggu     always process depths elementwise
c 18.05.1998	ggu	don't ask for description anymore
c 17.10.2001    ggu     accept also grd files with some missing data
c 18.10.2005    ggu     some error messages slightly changed
c 06.04.2009    ggu     read param.h
c 24.04.2009    ggu     new call to rdgrd()
c
c notes :
c
c could eliminate scrambling of iknot ->
c       no knscr()
c       pass ngrdim to bandop
c       change cmv,rosen
c
c********************************************************

        program vp

        include 'param.h'

	character*80 name
	character*80 file
	character*80 errfil
	character*80 errtex
        character*80 descrg,descrr,descra
        logical bstop,bopti

	common /descrr/descrr
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /pkonst/ grav,fcor,dcor,dirn

	common /xgv/xgv(nkndim), /ygv/ygv(nkndim)
	common /ipev/ipev(neldim), /ipv/ipv(nkndim)
	common /nen3v/nen3v(3,neldim), /hm3v/hm3v(3,neldim)
	common /iarv/iarv(neldim)

        integer ipdex(nkndim), iedex(neldim)

        integer iphv(nkndim), kphv(nkndim)
        integer iphev(neldim), iaux(neldim)
	integer neaux(3,neldim)
        real raux(neldim)
        real hev(neldim)
        real hkv(nkndim)
        real rphv(nkndim)

        integer ng(nkndim),iknot(ngrdim,nkndim)
        integer kvert(2,nkndim)

c        common /iphv/iphv(nkndim), /kphv/kphv(nkndim)
c        common /iphev/iphev(neldim), /iaux/iaux(neldim)
c        common /ng/ng(nkndim), /iknot/iknot(ngrdim,nkndim)
c        common /kvert/kvert(2,nkndim) !aux vector for rosen
c        common /hev/hev(neldim)
c        common /hkv/hkv(nkndim) !for depths

        integer nkn,nknh,nel,nelh,nli,nco
	integer nlidim,nlndim

	data bstop /.false./
	data errfil /'errout.dat'/
c
        dcor=0.
        dirn=0.
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
c
c        xscal=1.
c        yscal=1.
c        hscal=1.
	itief=0		!0=read by element  1=read by node
c
c open error file
c
	ner=ifileo(ner,errfil,'form','new')
	if(ner.le.0) then
		write(6,*) 'Cannot open error file'
		stop 'error stop : vp'
	end if
c
c get name of basin
c
	call pardef(0)
	ianz=igetxt('Enter name of basin : ',name)
	if(ianz.le.0) stop
	call putfnm('basnam',name)
c
c read all ?
c
c	idepth=-1
c	do while(idepth.eq.-1)
c	  idepth=iantw(' Change only in depth data ?')
c	  if(idepth.eq.-1) write(nat,*) ' Read error. Try again!'
c	end do
c
c always process whole file

	idepth = 0
c
c read input file
c
        file=name(1:ichanm(name))//'.grd'
        write(6,*) ' ...reading file ',file(1:ichanm(file))

	nlidim = 0
	nlndim = 0
        call rdgrd(
     +                   file
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iaux
     +                  ,iaux,iarv,iaux
     +                  ,hkv,hev,raux
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,iaux,iaux
     +                  )

c        call rdgrd(file,ner,bstop,nco,nkn,nknh,nel,nelh,nli
c     +                  ,nkndim,neldim,nliread
c     +                  ,ipv,ipev,iaux,iarv,nen3v,xgv,ygv,hev,hkv)

        write(nat,*) ' ...end of read'

	nknh = 0
	do k=1,nkn
	  if( hkv(k) .ne. -999. ) nknh = nknh + 1
	end do

	nelh = 0
	do k=1,nel
	  if( hev(k) .ne. -999. ) nelh = nelh + 1
	end do

        write(6,*) 'nkn,nel   : ',nkn,nel
        write(6,*) 'nknh,nelh : ',nknh,nelh
        write(6,*) 'nli,nco   : ',nli,nco

	if(bstop) then
	  goto 99999
        end if

        if(nkn.le.0 .or. nel.le.0) then
          write(ner,*) ' Nothing to process'
	  goto 99999
        end if

	if( nel.eq.nelh .and. nkn.eq.nknh ) then
	  write(nat,*) ' Can process depth node or elementwise.'
	  itief=0
          write(nat,*) ' ...depths are processed elementwise'
	else if(nel.eq.nelh) then
	  itief=0
          write(nat,*) ' ...depths are processed elementwise'
	else if(nkn.eq.nknh) then
	  itief=1
          write(nat,*) ' ...depths are processed nodewise'
	else if(nknh.eq.0.and.nelh.eq.0) then
	  itief=0
	  write(nat,*) ' No depth data read. Process anyway'
	else
	  itief=0
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	  write(nat,*) ' Not enough depth data read. Process anyway'
	  write(nat,*) '********************************************'
	  write(nat,*) '********************************************'
	end if
c
c process geometry ?
c
	if(idepth.ne.0) then
	  do k=1,nkn
	    iphv(k)=k
	  end do
	  do ie=1,nel
	    iphev(ie)=ie
	  end do
	  goto 1
	end if
c
c open files
c
        nb2=ideffi('basdir','basnam','.bas','unform','new')
        if(nb2.le.0) stop

c end reading ----------------------------------------------------

	call gtest('end read',neldim,nkn,nel,nen3v)

        bstop=.false.

        write(nat,*) ' ...making external index'

        call isort(nkn,ipv,ipdex)
        call isort(nel,ipev,iedex)

        write(nat,*) ' ...controlling uniqueness of node numbers'

        call uniqn(nkn,ipv,ipdex,ner,bstop)
	if(bstop) goto 99909

        write(nat,*) ' ...controlling uniqueness of element numbers'

        call uniqn(nel,ipev,iedex,ner,bstop)

	call gtest('end uniqn',neldim,nkn,nel,nen3v)

        write(nat,*) ' ...controlling uniqueness of elements'

        call uniqe(nel,nen3v,iaux,iphev,ipev,ner,bstop)
        if(bstop) goto 99915

	write(nat,*) ' ...changing extern with intern node numbers'

        call chexin(nkn,nel,nen3v,ipv,ipdex,ner,bstop)
	if(bstop) goto 99920

        write(nat,*) ' ...controlling node numbers'

        call needn(nkn,nel,nen3v,ipv,iaux,ner,bstop)
	if(bstop) goto 99918

	write(nat,*) ' ...testing sense of nodes in index'

	call gtest('end sense',neldim,nkn,nel,nen3v)

        call clockw(nkn,nel,nen3v,ipev,xgv,ygv,ner,bstop)
	if(bstop) goto 99930

	write(nat,*) ' ...setting up side index'

        call sidei(nkn,nel,nen3v,ng,iknot,ngrdim,ngr)
	call knscr(nkn,ngr,ngrdim,iknot)

	write(nat,*) ' Maximum grade of nodes is ',ngr

c bandwidth optimization -------------------------------------------

	call gtest('bandwidth',neldim,nkn,nel,nen3v)

        write(nat,*) ' ...optimizing band width'

        call bandw(nel,nen3v,mbw)

	call gtest('bandwidth 1',neldim,nkn,nel,nen3v)
	write(nat,*) ' Bandwidth is ',mbw

        call bandop(nkn,ngr,ipv,iphv,kphv,ng,iknot,kvert,bopti)

	call gtest('bandwidth 2',neldim,nkn,nel,nen3v)
        if(bopti) then
          call bandex(nkn,nel,nen3v,neaux,kphv,iphv,rphv
     +				,ipv,xgv,ygv,hkv)
	call gtest('bandwidth 3',neldim,nkn,nel,nen3v)
          call bandw(nel,nen3v,mbw)
	call gtest('bandwidth 4',neldim,nkn,nel,nen3v)
	  write(nat,*) ' Optimized bandwidth is ',mbw
	end if

c ------------------------------------------------------------------

	call gtest('end bandwidth',neldim,nkn,nel,nen3v)

	write(nat,*) ' ...renumbering elements'

        call renel(nel,nen3v,iaux,iedex,neaux,ipev,iarv,hev,raux)

c save pointers for depth ------------------------------------------

        write(nat,*) ' ...saving pointers'

	do k=1,nkn
	  iphv(k)=ipv(k)
	end do

	do ie=1,nel
	  iphev(ie)=ipev(ie)
	end do

c write to nb2 -----------------------------------------------------

	call gtest('write',neldim,nkn,nel,nen3v)

	call sp13uw(nb2)

	close(nb2)

c*****************************************
c	end part 1
c*****************************************

    1   continue
c
c open bas file
c
	nb2=ideffi('basdir','basnam','.bas','unform','old')
	if(nb2.le.0) stop
c
c read from nb2
c
	call sp13rr(nb2,nkndim,neldim)
c
	if(nkn.gt.nkndim) goto 99900
	if(nel.gt.neldim) goto 99900
	if(ngr.gt.ngrdim) goto 99900
c
	descrg=descrr
c
c process depths -------------------------------------------------

        write(nat,*) ' ...processing depths'
c
	if(itief.eq.0) then
          write(nat,*) ' ..........................elementwise'
          call helem(nel,nelh,iphev,iaux,iedex,ipev
     +				,hm3v,hev,ner,bstop)
	else
          write(nat,*) ' ..........................nodewise'
          call hnode(nkn,nel,nknh,nen3v,iphv,iaux,ipdex,ipv
     +                      ,hm3v,hkv,ner,bstop)
	end if

c
c get description of data file
c
c	write(nat,*)
c	write(nat,*) ' Possible descriptions for file :'
c	write(nat,*)
c        write(nat,*) '        1 --> (default)'
c	write(nat,*) descrg
c        write(nat,*) '        2 --> (new text)'
c	write(nat,*)
c
c        ianz=ideflt(1,' Enter value:')
c        if(ianz.eq.2) then
c		write(nat,*) ' Enter text :'
c		read(net,6000) descra
c		call tablnc(descra,descrr)
c	else
c          descrr=descrg
c	end if

c	write(nat,*) ' Description for file :'
c	write(nat,*) descrg

	descrr=descrg
c
c write to nb2
c
	write(nat,*) ' ...writing file'
	write(nat,*)
c
	call sp13uw(nb2)
c
	close(nb2)
c
	write(nat,*)
	write(nat,*) descrr
	write(nat,*)
	write(nat,*) ' nkn  = ',nkn ,'   nel  = ',nel
	write(nat,*) ' ngr  = ',ngr ,'   mbw  = ',mbw
	write(nat,*)
	write(nat,*) ' dcor = ',dcor,'   dirn = ',dirn
	write(nat,*)

	stop ' Successful completion'

99900	write(nat,*)' (00) error in dimension declaration'
	goto 99
c99901	write(nat,*)' (01) error reading first record of file 1'
c	goto 99
c99902	write(nat,*)' (02) error reading xscal,yscal'
c	goto 99
c99903	write(nat,*)' (03) error in version number :',nvers
c	goto 99
c99905   write(nat,*)' (05) error reading nanf,nend'
c        goto 99
c99907	write(nat,*)' (07) error reading coordinates close to'
c	write(nat,*)' node ',in
c	goto 99
99909	write(nat,*)' (09) error : no unique definition of nodes'
	goto 99
99915	write(nat,*)' (15) error : no unique definition of elements'
	goto 99
c99917	write(nat,*)' (17) error reading element index close to'
c	write(nat,*)' element ',nelm
c	goto 99
99918	write(nat,*)' (18) error reading file 1'
	goto 99
99920	write(nat,*)' (20) error reading file 1'
	goto 99
99930	write(nat,*)' (30) error in element index'
	goto 99
c99961	write(nat,*)' (61) error in version number :',nvers
c	goto 99
c99962	write(nat,*)' (62) error reading first record of file 3'
c	goto 99
c99963	write(nat,*)' (63) error reading hscal,itief'
c	goto 99
c99971	write(nat,*)' (71) error reading file 3'
c	goto 99
c99972   write(nat,*)' (72) error reading file 3'
c        goto 99
c99973	write(nat,*)' (73) error reading depth data close to'
c	write(nat,*)' node/element ',kn
c	goto 99
c99976   write(nat,*)' (76) error reading file 3'
c        goto 99
c99977	write(nat,*)' (77) error reading file 3'
c	goto 99
c99979   write(nat,*)' (79) error reading file 3'
c        goto 99
c99981	write(nat,*)' (81) error reading area code'
c	goto 99
c99982	write(nat,*)' (82) error reading area code close to'
c	write(nat,*)' nanf,nend,narea',nanf,nend,narea
c	goto 99
99999	write(nat,*)' (99) generic error'
	goto 99
c
   99	continue
c
	rewind(ner)
   98	read(ner,6000,err=97,end=97) errtex
	write(nat,*) errtex
	goto 98
   97	continue
c
	write(nat,*)
	write(nat,*) ' Error messages have been written to file :'
	write(nat,*) errfil
	write(nat,*)
	stop ' error stop : vpn2'
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
          do while(ie2.le.nel.and.iaux(index(ie1)).eq.iaux(index(ie2)))
            ie2=ie2+1
          end do
          if(ie2-ie1.gt.1) then
            isum=isum+ie2-ie1-1
            call equale(nel,ie1,ie2-1,nen3v,index,ipev,ner,bstop)
          end if
          ie1=ie2
        end do

        write(6,*) 'uniqe (isum) : ',isum

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

        subroutine sidei(nkn,nel,nen3v,ng,iknot,ngrdim,ngr)

c set up side index and find grade

        implicit none

        integer nkn,nel,ngr
        integer ngrdim
        integer nen3v(3,nel)
        integer ng(nkn)
        integer iknot(ngrdim,nkn)

        integer i,ii,iii,ie,k1,k2

	do i=1,nkn
	  do ii=1,ngrdim
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
	    if(ng(k1).gt.ngrdim) goto 99
            iknot(ng(k1),k1)=k2
            ng(k2)=ng(k2)+1        !has to be added to index
	    if(ng(k2).gt.ngrdim) goto 99
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
	write(6,*) 'dimension of ngrdim is too low: ',ngrdim
	stop 'error stop sidei: ngrdim'
        end

c**********************************************************

	subroutine knscr(nkn,ngr,ngrdim,iknot)

c scrambles side index for cmv and rosen

c -> please eliminate need for this

	implicit none

	integer nkn,ngr,ngrdim
        integer iknot(ngrdim*nkn)

	integer k,j

        do k=1,nkn
          do j=1,ngr
            iknot((k-1)*ngr+j)=iknot((k-1)*ngrdim+j)
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

        subroutine bandop(nkn,ngrdim,ipv,iphv,kphv,ng,iknot,kvert,bopti)

c optimize band width

        implicit none

        integer nkn,ngrdim
        logical bopti
        integer ipv(nkn)
        integer iphv(nkn),kphv(nkn)
        integer ng(nkn)
        integer iknot(ngrdim,nkn)
        integer kvert(2,nkn)

        integer iantw
	logical bauto

        bopti = iantw(' Optimization of bandwidth ?') .gt. 0
	if( .not. bopti ) return

        bauto = iantw(' Automatic optimization ?') .gt. 0
	if( bauto ) then
	  call ininum(nkn,iphv,kphv)
	  call optest(nkn,iphv,kphv)
	  call cmgrade(nkn,ngrdim,ipv,iphv,kphv,ng,iknot,1,4)
	  call optest(nkn,iphv,kphv)
	  call revnum(nkn,iphv,kphv)
	  call optest(nkn,iphv,kphv)
          call rosen(nkn,ngrdim,iphv,kphv,ng,iknot,kvert)
	  call optest(nkn,iphv,kphv)
	  return
	end if

        do while( bopti )
	  call ininum(nkn,iphv,kphv)

c	  call anneal(nkn,ngrdim,kphv,ng,iknot,iphv,kvert)

          if( iantw(' Cuthill McKee algorithm ?') .gt. 0 ) then
            call cmv(nkn,ngrdim,ipv,iphv,kphv,ng,iknot)
	  end if

          if( iantw(' Reverse numbering of nodes ?') .gt. 0 ) then
	    call revnum(nkn,iphv,kphv)
	  end if

          if( iantw(' Rosen algorithm ?') .gt. 0 ) then
            call rosen(nkn,ngrdim,iphv,kphv,ng,iknot,kvert)
	  end if

	  bopti = iantw(' Repeat optimization of bandwidth ?') .gt. 0
        end do

	bopti = .true.

        end

c**********************************************************

	subroutine revnum(nkn,iphv,kphv)

c reverses numbering of nodes

	implicit none

	integer nkn
	integer iphv(1), kphv(1)

	integer i

	write(6,*) 'Applying reverse algorithm...'

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
	integer iphv(1), kphv(1)

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
	integer iphv(1), kphv(1)

	integer i

	do i=1,nkn
	  iphv(i) = 0
	  kphv(i) = 0
	end do

	end

c**********************************************************

        subroutine bandex(nkn,nel,nen3v,neaux,kphv,iphv,rphv
     +                      ,ipv,xgv,ygv,hkv)

c exchange nodes after optimization

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer neaux(3,nel)
        integer kphv(nkn)
        integer ipv(nkn)
        integer iphv(nkn)
        real rphv(nkn)
        real xgv(nkn),ygv(nkn)
        real hkv(nkn)

        integer ie,ii

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

c       copy arrays with kphv as rank table (iphv,rphv are aux arrays)

        call icopy(nkn,ipv,iphv,kphv)
        call rcopy(nkn,xgv,rphv,kphv)
        call rcopy(nkn,ygv,rphv,kphv)
        call rcopy(nkn,hkv,rphv,kphv)

        return
        end

c**********************************************************

        subroutine icopy(n,iv,iauxv,irank)

c copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        integer iv(n),iauxv(n)
        integer irank(n)

        integer i

        do i=1,n
          iauxv(i)=iv(i)
        end do

        do i=1,n
          iv(irank(i))=iauxv(i)
        end do

        return
        end

c**********************************************************

        subroutine rcopy(n,rv,rauxv,irank)

c copy one array to itself exchanging elements as in irank

        implicit none

        integer n
        real rv(n),rauxv(n)
        integer irank(n)

        integer i

        do i=1,n
          rauxv(i)=rv(i)
        end do

        do i=1,n
          rv(irank(i))=rauxv(i)
        end do

        return
        end

c**********************************************************

        subroutine renel(nel,nen3v,iaux,iedex,neaux,ipev,iarv,hev,raux)

c renumbering of elements

c we construct iedex newly and use iaux,iedex as aux arrays
c iedex is also used as a real aux array for rcopy	-> changed
c neaux is probably he3v (real) used as an aux array	-> changed

        implicit none

        integer nel
        integer nen3v(3,nel)
        integer iaux(nel)
        integer iedex(nel)
        integer neaux(3,nel)
        integer ipev(nel)
	integer iarv(nel)
        real hev(nel)
        real raux(nel)

        integer ie,ii

        do ie=1,nel
          iaux(ie)=min(nen3v(1,ie),nen3v(2,ie),nen3v(3,ie))
	end do

        call isort(nel,iaux,iedex)  !iedex is the index table
        call rank(nel,iedex,iaux)   !iaux is the rank table

c       now we use iedex as an aux array and sort with iaux (rank)

        call icopy(nel,ipev,iedex,iaux)
        call icopy(nel,iarv,iedex,iaux)
        call rcopy(nel,hev,raux,iaux)  !we use iedex as real aux array

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,iaux(ie))=neaux(ii,ie)
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

	write(6,*) 'helem: ',nel,nhd

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
            do ii=1,3
              hm3v(ii,ie)=-999.
            end do
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
            if(iaux(kn).gt.0) then
              hm3v(ii,ie)=hkv(iaux(kn))
	    else
              hm3v(ii,ie)=-999.
            end if
          end do
	end do

        return
        end

c**********************************************************

	subroutine gtest(text,neldim,nkn,nel,nen3v)

	implicit none

	character*(*) text
	integer neldim,nkn,nel
	integer nen3v(3,neldim)

	integer ie,ii,k
	integer iii
	logical berror
	
	berror = .false.

	write(6,*) 'testing ... ',text

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

	subroutine optest(nkn,iphv,kphv)

	implicit none

	integer nkn
	integer iphv(1)
	integer kphv(1)

	integer i

	do i=1,nkn
	  if( iphv(i) .gt. nkn ) goto 99
	  if( kphv(i) .gt. nkn ) goto 99
	  if( iphv(kphv(i)) .ne. i ) goto 99
	  if( kphv(iphv(i)) .ne. i ) goto 99
	end do

	return
   99	continue
	write(6,*) 'error in pointers...'
	write(6,*) i,nkn
	write(6,*) iphv(i),kphv(i)
	stop 'error stop optest: pointers'
	end

c**********************************************************


