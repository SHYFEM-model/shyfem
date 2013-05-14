c
c $Id: subboxa.f,v 1.25 2009-05-21 09:24:00 georg Exp $
c
c subroutines for computing discharge / flux for boxes
c
c contents :
c
c subroutine inboxa
c subroutine rdboxa
c subroutine ckboxa
c subroutine prboxa
c subroutine tsboxa
c
c subroutine wrboxa(it)				write of flux data
c
c subroutine boxini				initializes flux routines
c
c revision log :
c
c 30.04.1998    ggu	newly written routines (subpor deleted)
c 07.05.1998    ggu	check nrdveci on return for error
c 08.05.1998    ggu	restructured with new comodity routines
c 13.09.1999    ggu	type of node computed in own routine flxtype
c 19.11.1999    ggu	iskadj into sublin
c 20.01.2000	ggu	old routines substituted, new routine extrsect
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 20.01.2000    ggu     common block /dimdim/ eliminated
c 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
c 26.05.2003	ggu	in flxnov substituted a,b with b,c
c 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
c 10.08.2003	ggu	do not call setweg, setnod, setkan
c 23.03.2006    ggu     changed time step to real
c 28.09.2007    ggu     use testbndo to determine boundary node in flxtype
c 28.04.2009    ggu     links re-structured
c 23.02.2011    ggu     new routine call write_node_fluxes() for special output
c 01.06.2011    ggu     documentation to flxscs() changed
c 21.09.2011    ggu     some lower-level subroutines copied to subflx.f
c 07.10.2011    ggu     adjusted for 3d flux routines
c 19.10.2011    ggu     added T/S variables, created fluxes_*() routines
c 19.10.2011    ggu     added conz variables, created fluxes_template()
c 10.05.2013    ggu     adapted to boxes computations
c 14.05.2013    ggu     write also OBC sections and contributions
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

        subroutine mod_box(mode)
 
        implicit none
 
        integer mode
 
        include 'modules.h'
 
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
 
        if( mode .eq. M_AFTER ) then
           call wrboxa(it)
        else if( mode .eq. M_INIT ) then
           call inboxa
        else if( mode .eq. M_READ ) then
           call rdboxa
        else if( mode .eq. M_CHECK ) then
           call ckboxa
        else if( mode .eq. M_SETUP ) then
           call boxini
        else if( mode .eq. M_PRINT ) then
           call prboxa
        else if( mode .eq. M_TEST ) then
           call tsboxa
        else if( mode .eq. M_BEFOR ) then
c          nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_box'
        end if
 
        end

c******************************************************************

        subroutine inboxa

c nsect		total number of sections
c kfluxm	total number of nodes defining sections
c kflux()	node numbers defining sections

        implicit none

	include 'param.h'
	include 'subboxa.h'

        nsect = -1	!must still be initialized
        kfluxm = 0

        end

c******************************************************************

        subroutine rdboxa

        implicit none

	include 'param.h'
	include 'subboxa.h'

	integer ndim
        integer nrdveci

	ndim = nfxboxdim

        kfluxm = nrdveci(kflux,ndim)

        if( kfluxm .lt. 0 ) then
          if( kfluxm .eq. -1 ) then
            write(6,*) 'dimension error nfxdim in section $flux : '
     +                          ,ndim
          else
            write(6,*) 'read error in section $fbox'
          end if
          stop 'error stop rdboxa'
        end if

        end

c******************************************************************

        subroutine ckboxa

        implicit none

	include 'param.h'
	include 'subboxa.h'

	integer k,ii
        logical berror

c convert external to internal node numbers

	call n2int(kfluxm,kflux,berror)

        if( berror ) then
		write(6,*) 'error in section FLUX'
		stop 'error stop: ckboxa'
	end if

c the real set up is done in boxini
c but since at this stage we do not have all the arrays set up
c we post-pone it until later

        end

c******************************************************************

	subroutine prboxa

	implicit none

	include 'param.h'
	include 'subboxa.h'
	
	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii

	integer ipext
	logical nextline

	write(6,*)
	write(6,*) 'flux section :'
	write(6,*)
	write(6,*) 'nsect,kfluxm ',nsect,kfluxm
	write(6,*)

	ns = 0
	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'section : ',ns,ntotal
	  do i=ifirst,ilast
	    write(6,*) ipext(kflux(i)),(iflux(ii,i),ii=1,3)
	  end do
	end do

	end

c******************************************************************

	subroutine tsboxa

	implicit none

	include 'param.h'
	include 'subboxa.h'
	
	integer i,ii

	write(6,*) '/kfluxc/'
	write(6,*) nsect,kfluxm
	write(6,*) (kflux(i),i=1,kfluxm)

	write(6,*) '/iflux/'
	write(6,*) ((iflux(ii,i),ii=1,3),i=1,kfluxm)

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine wrboxa(it)

c administers writing of flux data

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer it

	integer itanf,itend
	integer j,i,l,lmax,nlmax,ivar,nvers
	real az,azpar,rr

	real rhov(nlvdim,nkndim)
        common /rhov/rhov
        real saltv(nlvdim,nkndim)
        common /saltv/saltv
        real tempv(nlvdim,nkndim)
        common /tempv/tempv
        real cnv(nlvdim,nkndim)
        common /cnv/cnv
        real zenv(3,neldim)
        common /zenv/zenv

	integer ifemop
	real getpar

	real fluxes(0:nlvdim,3,nscboxdim)

	integer nlayers(nscboxdim)		!number of layers in section
	real masst(0:nlvdim,3,nscboxdim)	!accumulator mass
	real saltt(0:nlvdim,3,nscboxdim)	!accumulator salt
	real tempt(0:nlvdim,3,nscboxdim)	!accumulator temp
	real conzt(0:nlvdim,3,nscboxdim)	!accumulator conz

	integer nlayers_ob(nbcdim)		!number of layers in section
	real fluxes_ob(0:1,3,nbcdim)		!discharges at open bound
	real masst_ob(0:1,3,nbcdim)		!discharges at open bound

	integer nnt,nns,nne,nnob
	save nnt,nns,nne,nnob
	real valt(nbxdim)
	real vals(nbxdim)
	real vale(nbxdim)
	real val(nbxdim)
	save valt,vals,vale

	save nlayers
        save masst,saltt,tempt,conzt

        integer nrm,nrs,nrt,nrc
	save nrm,nrs,nrt,nrc

        integer idtbox,itbox,itmbox
        save idtbox,itbox,itmbox
        integer nbbox
        save nbbox
	integer ibarcl,iconz
	save ibarcl,iconz

	integer nbnds

        data nbbox /0/

c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------

        if( nbbox .eq. -1 ) return

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        if( nbbox .eq. 0 ) then

                idtbox = nint(getpar('idtbox'))
                itmbox = nint(getpar('itmbox'))
                itanf = nint(getpar('itanf'))
                itend = nint(getpar('itend'))

                if( idtbox .le. 0 ) nbbox = -1
                if( itmbox .gt. itend ) nbbox = -1
                if( nbbox .eq. -1 ) return

		call box_init

                if( kfluxm .le. 0 ) nbbox = -1
                if( nsect .le. 0 ) nbbox = -1
                if( nbbox .eq. -1 ) return

                if( nsect .gt. nscboxdim ) then
                  stop 'error stop wrboxa: dimension nscboxdim'
                end if

		ibarcl = nint(getpar('ibarcl'))
		iconz = nint(getpar('iconz'))
		ibarcl = 0
		iconz = 0

		if( itmbox .lt. itanf ) itmbox = itanf
                itbox = itmbox + idtbox
		itmbox = itmbox + 1	!start from next time step

		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		call fluxes_init(nlvdim,nsect,nlayers,nrm,masst)
		if( ibarcl .gt. 0 ) then
		  call fluxes_init(nlvdim,nsect,nlayers,nrs,saltt)
		  call fluxes_init(nlvdim,nsect,nlayers,nrt,tempt)
		end if
		if( iconz .eq. 1 ) then
		  call fluxes_init(nlvdim,nsect,nlayers,nrc,conzt)
		end if

		call boxes_init(nbox,nnt,valt)	!temp in box
		call boxes_init(nbox,nns,vals)	!salt in box
		call boxes_init(nbox,nne,vale)	!zeta in box

		nbc_ob = nbnds()
		do i=1,nbc_ob
		  nlayers_ob(i) = 1
		end do
		call fluxes_init(1,nbc_ob,nlayers_ob,nnob,masst_ob)

                nbbox=ifemop('.box','unform','new')
                if(nbbox.le.0) then
        	   stop 'error stop wrboxa : Cannot open FLX file'
		end if

	        nvers = 5
                call wfflx      (nbbox,nvers
     +                          ,nsect,kfluxm,idtbox,nlmax
     +                          ,kflux
     +                          ,nlayers
     +                          )

c               here we could also compute and write section in m**2

        end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

        if( it .lt. itmbox ) return

	call getaz(azpar)
	az = azpar

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)
	call fluxes_accum(nlvdim,nsect,nlayers,nrm,masst,fluxes)

c	write(6,*) (nlayers(i),i=1,nsect)
c	write(6,*) (((masst(l,j,i),l=0,nlayers(i)),j=1,3)
c     +                  ,i=1,nsect)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,saltv)
	  call fluxes_accum(nlvdim,nsect,nlayers,nrs,saltt,fluxes)
	  ivar = 12
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,tempv)
	  call fluxes_accum(nlvdim,nsect,nlayers,nrt,tempt,fluxes)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,cnv)
	  call fluxes_accum(nlvdim,nsect,nlayers,nrc,conzt,fluxes)
	end if

	call box_scalar(tempv,val)
	call boxes_accum(nbox,nnt,valt,val)
	call box_scalar(saltv,val)
	call boxes_accum(nbox,nns,vals,val)
	call box_eta(val)
	call boxes_accum(nbox,nne,vale,val)

	call box_ob_compute(fluxes_ob)	!open boundary conditions
	call fluxes_accum(1,nbc_ob,nlayers_ob,nnob,masst_ob,fluxes_ob)

c	-------------------------------------------------------
c	time for output?
c	-------------------------------------------------------

        if( it .lt. itbox ) return
        itbox=itbox+idtbox

c	-------------------------------------------------------
c	average and write results
c	-------------------------------------------------------

	ivar = 0
	call fluxes_aver(nlvdim,nsect,nlayers,nrm,masst,fluxes)
	call wrflx(nbbox,it,nlvdim,nsect,ivar,nlayers,fluxes)

	if( ibarcl .gt. 0 ) then
	  ivar = 11
	  call fluxes_aver(nlvdim,nsect,nlayers,nrs,saltt,fluxes)
	  call wrflx(nbbox,it,nlvdim,nsect,ivar,nlayers,fluxes)
	  ivar = 12
	  call fluxes_aver(nlvdim,nsect,nlayers,nrt,tempt,fluxes)
	  call wrflx(nbbox,it,nlvdim,nsect,ivar,nlayers,fluxes)
	end if

	if( iconz .eq. 1 ) then
	  ivar = 10
	  call fluxes_aver(nlvdim,nsect,nlayers,nrc,conzt,fluxes)
	  call wrflx(nbbox,it,nlvdim,nsect,ivar,nlayers,fluxes)
	end if

	call boxes_aver(nbox,nnt,valt)	!temp in box
	call boxes_aver(nbox,nns,vals)	!salt in box
	call boxes_aver(nbox,nne,vale)	!zeta in box

	call fluxes_aver(1,nbc_ob,nlayers_ob,nnob,masst_ob,fluxes_ob)

	call box_write(it,fluxes,valt,vals,vale,fluxes_ob)

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

	call fluxes_init(nlvdim,nsect,nlayers,nrm,masst)

	if( ibarcl .gt. 0 ) then
	  call fluxes_init(nlvdim,nsect,nlayers,nrs,saltt)
	  call fluxes_init(nlvdim,nsect,nlayers,nrt,tempt)
	end if

	if( iconz .eq. 1 ) then
	  call fluxes_init(nlvdim,nsect,nlayers,nrc,conzt)
	end if

	call boxes_init(nbox,nnt,valt)	!temp in box
	call boxes_init(nbox,nns,vals)	!salt in box
	call boxes_init(nbox,nne,vale)	!zeta in box

	call fluxes_init(1,nbc_ob,nlayers_ob,nnob,masst_ob)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine boxini

c initializes flux routines finally (wrapper for flx_init)

	implicit none

	include 'param.h'
	include 'subboxa.h'
	
	call flx_init(kfluxm,kflux,nsect,iflux)

	end

c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_init

c reads file etc...

	implicit none

	call inboxa		!basic initialization

	call box_read		!reads boxes.txt
	call box_elab		!sets up ikboxes

	call ckboxa		!converts ext to int nodes

	call box_ob_init	!initializes OB contributions (nodes intern)

	call boxini		!sets data structure iflux

	end

c******************************************************************

	subroutine box_read

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'subboxa.h'

	integer nb,id,n,ib1,ib2,i
	integer nelaux,ie
	integer ierr,iaux
	integer iauxv(nkndim)

	nsect = 0
	kfluxm = 0
	ierr = 0

	write(6,*) 'start reading boxes.txt... '
	open(1,file='boxes.txt',form='formatted',status='old',err=94)

	read(1,*) nelaux,nbox,nb
	if( nelaux .ne. nel ) goto 99
	if( nbox .gt. nbxdim ) goto 98
	read(1,*) (iboxes(ie),ie=1,nelaux)

    2	continue
	  read(1,*) id,n,ib1,ib2
	  if( id .le. 0 ) goto 3
	  read(1,*) (iauxv(i),i=1,n)
	  nsect = nsect + 1
	  if( id .ne. nsect ) goto 97

	  call box_insert_section(id,n,ib1,ib2,iauxv,ierr)

	  goto 2
    3	continue

	close(1)

	if( nsect .gt. nscboxdim ) goto 95
	if( kfluxm .gt. nfxboxdim ) goto 96

	kfluxm = kfluxm - 1

	write(6,*) 'finished reading boxes.txt: '
	write(6,*) nbox,nsect,kfluxm

	return
   94	continue
	write(6,*) 'file: boxes.txt'
	stop 'error stop box_read: cannot read file'
   95	continue
	write(6,*) nsect,nscboxdim
	stop 'error stop box_read: nscboxdim'
   96	continue
	write(6,*) kfluxm,nfxboxdim
	stop 'error stop box_read: nfxboxdim'
   97	continue
	write(6,*) id,nsect
	stop 'error stop box_read: id'
   98	continue
	write(6,*) nbox,nbxdim
	stop 'error stop box_read: nbxdim'
   99	continue
	write(6,*) nel,nelaux
	stop 'error stop box_read: nel'
	end

c******************************************************************

	subroutine box_insert_section(id,n,ib1,ib2,list,ierr)

c inserts section in data structure

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'subboxa.h'

	integer id,n,ib1,ib2
	integer list(n)
	integer ierr

	integer i

	if( id .gt. nscboxdim ) ierr = ierr + 1
	if( kfluxm + n + 1 .gt. nfxboxdim ) ierr = ierr + 1

	if( ierr .eq. 0 ) then	! no errors -> insert
	    isects(1,id) = n
	    isects(2,id) = kfluxm+1
	    isects(3,id) = ib1
	    isects(4,id) = ib2
	    do i=1,n
	      kflux(kfluxm+i) = list(i)
	    end do
	    kfluxm = kfluxm + n + 1
	    kflux(kfluxm) = 0
	else
	    kfluxm = kfluxm + n + 1
	end if

	end

c******************************************************************

	subroutine box_elab

c sets up ikboxes - index of nodes to boxes
c
c is -1 if node belongs to more than one box

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'subboxa.h'

	integer k,ie,ii,ib
	integer nb,nl

	do k=1,nkn
	  ikboxes(k) = 0
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ikboxes(k) .eq. 0 ) then
	      ikboxes(k) = ib
	    else if( ikboxes(k) .ne. ib ) then
	      ikboxes(k) = -1
	    end if
	  end do
	end do

	nb = 0
	nl = 0
	do k=1,nkn
	  if( ikboxes(k) .gt. 0 ) then
	    nb = nb + 1
	  else if( ikboxes(k) .lt. 0 ) then
	    nl = nl + 1
	  end if
	end do

	if( nb+nl .ne. nkn ) stop 'error stop box_elab: internal error'

	write(6,*) 'box nodes: ',nb,nl,nb+nl,nkn

	end

c******************************************************************

	subroutine box_write(it,fluxes,valt,vals,vale,masst_ob)

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'subboxa.h'

	integer it
	real fluxes(0:nlvdim,3,nscboxdim)
	real valt(nbxdim)
	real vals(nbxdim)
	real vale(nbxdim)
	real masst_ob(0:1,3,nbcdim)

	integer iu,ib,is,ib1,ib2,ii

	iu = 77

	write(iu,*) it

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,valt(ib),vals(ib),vale(ib)
	end do

	write(iu,*) nsect
	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  write(iu,*) ib1,ib2,(fluxes(0,ii,is),ii=1,3)
	end do

	write(iu,*) nbc_ob
	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  write(iu,*) ib1,ib2,(masst_ob(0,ii,is),ii=1,3)
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_eta(val)

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'ev.h'
	include 'subboxa.h'

	real val(nbxdim)

        real zenv(3,neldim)
        common /zenv/zenv
        integer iwegv(neldim)
        common /iwegv/iwegv

	double precision vald(nbxdim)
	double precision vold(nbxdim)

	integer ib,ie,ii
	real vol,area3,z

	do ib=1,nbox
	  vald(ib) = 0.
	  vold(ib) = 0.
	end do

	do ie=1,nel
	  if( iwegv(ie) .eq. 0 ) then	!only for wet elements
	    ib = iboxes(ie)
	    area3 = 4.*ev(10,ie)
	    do ii=1,3
	      z = zenv(ii,ie)
	      vol = area3
	      vald(ib) = vald(ib) + vol*z
	      vold(ib) = vold(ib) + vol
	    end do
	  end if
	end do

	do ib=1,nbox
	  if( vold(ib) .gt. 0. ) then
	    val(ib) = vald(ib)/vold(ib)
	  else
	    val(ib) = 0.
	  end if
	end do

	end

c******************************************************************

	subroutine box_scalar(scalar,val)

	implicit none

	include 'param.h'
	include 'basin.h'
	include 'ev.h'
	include 'subboxa.h'

	real scalar(nlvdim,nkndim)
	real val(*)

        integer ilhv(neldim)
        common /ilhv/ilhv
        real hdknv(nlvdim,nkndim)
        common /hdknv/hdknv

	double precision vald(nbxdim)
	double precision vold(nbxdim)

	integer ib,ie,lmax,l,ii,k
	real vol,area3

	do ib=1,nbox
	  vald(ib) = 0.
	  vold(ib) = 0.
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vald(ib) = vald(ib) + vol*scalar(l,k)
	      vold(ib) = vold(ib) + vol
	    end do
	  end do
	end do

	do ib=1,nbox
	  if( vold(ib) .gt. 0. ) then
	    val(ib) = vald(ib)/vold(ib)
	  else
	    val(ib) = 0.
	  end if
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine boxes_init(n,na,val)

	implicit none

	integer n,na
	real val(n)

	integer i

	na = 0
	do i=1,n
	  val(i) = 0
	end do

	end

c******************************************************************

	subroutine boxes_accum(n,na,val,value)

	implicit none

	integer n,na
	real val(n)
	real value(n)

	integer i

	na = na + 1
	do i=1,n
	  val(i) = val(i) + value(i)
	end do

	end

c******************************************************************

	subroutine boxes_aver(n,na,val)

	implicit none

	integer n,na
	real val(n)

	integer i

	if( na .le. 0 ) return

	do i=1,n
	  val(i) = val(i) / na
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_ob_init

c sets up open boundary inputs

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer nbc,ibc,itype
	integer nk,i,k,ib,ibk
	integer ierr
	integer iauxv(nkndim)

	integer nbnds,itybnd,nkbnds,kbnds
	integer ipext

	ierr = 0
	kfluxm = kfluxm + 1
	kflux(kfluxm) = 0

	nbc = nbnds()

	do ibc = 1,nbc
	  itype = itybnd(ibc)
	  nk = nkbnds(ibc)
	  ib = 0
	  do i=1,nk
	    k = kbnds(ibc,i)
	    ibk = ikboxes(k)
	    if( ibk .le. 0 ) goto 99	!node not in unique box
	    if( ib .eq. 0 ) ib = ibk
	    if( ibk .ne. ib ) goto 98	!boundary nodes in different boxes
	  end do
	  iscbnd(1,ibc) = nk		!total number of boundary nodes
	  iscbnd(2,ibc) = itype		!type of boundary
	  iscbnd(3,ibc) = -ibc		!from box
	  iscbnd(4,ibc) = ib		!to box

	  if( itype .eq. 1 ) then	!insert sections of z boundary
	    nsect = nsect + 1
	    call irbnds(ibc,nkndim,nk,iauxv)
	    call box_insert_section(nsect,nk,-ibc,ib,iauxv,ierr)
	  end if
	end do

	if( nsect .gt. nscboxdim ) goto 95
	if( kfluxm .gt. nfxboxdim ) goto 96

	kfluxm = kfluxm - 1

	return
   95	continue
	write(6,*) nsect,nscboxdim
	stop 'error stop box_ob_init: nscboxdim'
   96	continue
	write(6,*) kfluxm,nfxboxdim
	stop 'error stop box_ob_init: nfxboxdim'
   98	continue
	write(6,*) 'boxes: ',ib,ibk,'  node: ',ipext(k)
	stop 'error stop box_ob_init: different boxes for nodes'
   99	continue
	write(6,*) 'boxes: ',ibk,'  node: ',ipext(k)
	stop 'error stop box_ob_init: not unique box for node'
	end

c******************************************************************

	subroutine box_ob_compute(fluxes_ob)

c computes open boundary inputs

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real fluxes_ob(0:1,3,nbcdim)		!discharges at open bound

	integer ibc,nbc
	real val

	real get_discharge

	nbc = nbc_ob

	do ibc = 1,nbc
	  val = get_discharge(ibc)
	  if( val .ge. 0 ) then
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = val
	    fluxes_ob(0,3,ibc) = 0.
	    fluxes_ob(1,1,ibc) = val
	    fluxes_ob(1,2,ibc) = val
	    fluxes_ob(1,3,ibc) = 0.
	  else
	    fluxes_ob(0,1,ibc) = val
	    fluxes_ob(0,2,ibc) = 0.
	    fluxes_ob(0,3,ibc) = -val
	    fluxes_ob(1,1,ibc) = val
	    fluxes_ob(1,2,ibc) = 0.
	    fluxes_ob(1,3,ibc) = -val
	  end if
	end do

	end

c******************************************************************

