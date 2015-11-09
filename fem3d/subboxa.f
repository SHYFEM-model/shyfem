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
c 29.10.2014    ggu     new code and 3d version
c
c notes :
c
c insert parameters idtbox and itmbox into STR file to have box file written
c
c format of written files  please see box_write_*()
c
c still to check: some sections have more layers than adjacent boxes
c
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************
c******************************************************************

        subroutine mod_box(mode)
 
c this routine is not executed

        implicit none
 
        integer mode
 
        include 'modules.h'
 
	include 'femtime.h'
 
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
            write(6,*) 'dimension error ndim in section $flux : '
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

	use mod_conz, only : cnv
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv

	implicit none

	include 'param.h'
	include 'subboxa.h'
	!include 'femtime.h'

	integer it

	integer j,i,l,lmax,nlmax,ivar,nvers
	integer date,time
	integer idtbox
	integer itanf,itend
	real az,azpar,rr
	real dt



	integer ifemop
	real getpar
	double precision dgetpar
	logical has_output,is_in_output,next_output

	real fluxes(0:nlvdim,3,nscboxdim)
	real val2d(nbxdim)
	real val3d(0:nlvdim,nbxdim)

        integer nrm,nrs,nrt,nrc
	save nrm,nrs,nrt,nrc
	integer nlayers(nscboxdim)		!number of layers in section
	save nlayers
	real masst(0:nlvdim,3,nscboxdim)	!accumulator mass
	real saltt(0:nlvdim,3,nscboxdim)	!accumulator salt
	real tempt(0:nlvdim,3,nscboxdim)	!accumulator temp
	real conzt(0:nlvdim,3,nscboxdim)	!accumulator conz
        save masst,saltt,tempt,conzt

	integer nnob
	save nnob
	integer nlayers_ob(nbcdim)		!number of layers in section
	real fluxes_ob(0:1,3,nbcdim)		!discharges OBC (aux)
	real masst_ob(0:1,3,nbcdim)		!discharges OBC (accum)
	save nlayers_ob
        save fluxes_ob,masst_ob

	real dtbox
	save dtbox
	integer nnn
	save nnn

	integer nblayers(nbxdim)		!number of layers in boxes
	real barea(nbxdim)			!area of boxes
	real bvolume(nbxdim)			!volume of boxes
	real bdepth(nbxdim)			!depth of boxes
	save nblayers,barea,bvolume,bdepth

	real valt(0:nlvdim,nbxdim)
	real vals(0:nlvdim,nbxdim)
	real valv(0:nlvdim,nbxdim)
	real vale(nbxdim)
	real valm(7,nbxdim)
	save valt,vals,vale,valv,valm

	real eta_old(nbxdim)
	real eta_act(nbxdim)
	save eta_old,eta_act

	real valw(0:nlvdim,nbxdim)
	real valvis(0:nlvdim,nbxdim)
	real valdif(0:nlvdim,nbxdim)
	save valw,valvis,valdif

        integer nbbox
        save nbbox
	integer ibarcl,iconz,ievap
	save ibarcl,iconz,ievap
	integer ia_out(4)
	save ia_out

	integer nbnds
	logical bbox3d,bfirst
	save bfirst
	data bfirst /.true./

        data nbbox /0/

c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------

        if( nbbox .eq. -1 ) return

	bbox3d = .false.	!write 3d results if in 3d
	bbox3d = .true.		!write 3d results if in 3d

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        if( nbbox .eq. 0 ) then

          	call init_output('itmbox','idtbox',ia_out)
		call increase_output(ia_out)  !itbox=itmbox+idtbox
          	if( .not. has_output(ia_out) ) nbbox = -1

                if( nbbox .eq. -1 ) return

                date = nint(dgetpar('date'))
                time = nint(dgetpar('time'))

		call box_init
		call box_make_stats(nblayers,barea,bvolume,bdepth)

                if( kfluxm .le. 0 ) nbbox = -1
                if( nsect .le. 0 ) nbbox = -1
                if( nbbox .eq. -1 ) return

                if( nsect .gt. nscboxdim ) then
                  stop 'error stop wrboxa: dimension nscboxdim'
                end if

		ibarcl = nint(getpar('ibarcl'))
		iconz = nint(getpar('iconz'))
		ievap = nint(getpar('ievap'))
		ibarcl = 0
		iconz = 0

		call get_nlayers(kfluxm,kflux,nlayers,nlmax)

		call fluxes_init(nlvdim,nsect,nlayers,nrm,masst)
		if( ibarcl .gt. 0 ) then
		  call fluxes_init(nlvdim,nsect,nlayers,nrs,saltt)
		  call fluxes_init(nlvdim,nsect,nlayers,nrt,tempt)
		end if
		if( iconz .eq. 1 ) then
		  call fluxes_init(nlvdim,nsect,nlayers,nrc,conzt)
		end if

		nnn = 0
		dtbox = 0.
		call boxes_3d_init(nbox,valt)	!temp in box
		call boxes_3d_init(nbox,vals)	!salt in box
		call boxes_2d_init(nbox,vale)	!zeta in box
		call boxes_3d_init(nbox,valv)	!vel speed in box
		call boxes_multi_init(nbox,7,valm)	!meteo in box

		call boxes_3d_init(nbox,valw)	!vertical vel
		call boxes_3d_init(nbox,valvis)	!viscosity
		call boxes_3d_init(nbox,valdif)	!diffusivity

		call box_eta(zeov,val2d)
		eta_old = val2d

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
		idtbox = ia_out(1)
                call wfflx      (nbbox,nvers
     +                          ,nsect,kfluxm,idtbox,nlmax
     +                          ,kflux
     +                          ,nlayers
     +                          )

c               here we could also compute and write section in m**2

		itanf = nint(dgetpar('itanf'))
		itend = nint(dgetpar('itend'))
		call box_write_stats(itanf,itend,date,time
     +					,nlayers,nblayers
     +					,barea,bvolume,bdepth)

        end if

c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

        if( .not. is_in_output(ia_out) ) return		! it < itmbox

	if( bfirst ) then	!initial condition
	  bfirst = .false.
	  call box_eta(zenv,val2d)
	  eta_old = val2d
	  call box_write_init(eta_old)
	  return		!initial call to wrboxa
	end if

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,rhov)
	call fluxes_accum(nlvdim,nsect,nlayers,nrm,masst,fluxes)

	call box_ob_compute(fluxes_ob)	!open boundary conditions
	call fluxes_accum(1,nbc_ob,nlayers_ob,nnob,masst_ob,fluxes_ob)

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

	nnn = nnn + 1
	dtbox = dtbox + dt

	call box_3d_scalar(tempv,val3d)
	call boxes_3d_accum(nbox,valt,val3d)
	call box_3d_scalar(saltv,val3d)
	call boxes_3d_accum(nbox,vals,val3d)
	call box_vel(val3d)
	call boxes_3d_accum(nbox,valv,val3d)
	call box_eta(zenv,val2d)
	call boxes_2d_accum(nbox,vale,val2d)
	eta_act = val2d

	call box_3d_vertical(wlnv,val3d,0)
	call boxes_3d_accum(nbox,valw,val3d)
	call box_3d_vertical(visv,val3d,1)
	call boxes_3d_accum(nbox,valvis,val3d)
	call box_3d_vertical(difv,val3d,1)
	call boxes_3d_accum(nbox,valdif,val3d)

	call boxes_meteo_accum(nbox,valm,val2d,ievap)

c	-------------------------------------------------------
c	time for output?
c	-------------------------------------------------------

        if( .not. next_output(ia_out) ) return

c	-------------------------------------------------------
c	average results
c	-------------------------------------------------------

	ivar = 0
	call fluxes_aver(nlvdim,nsect,nlayers,nrm,masst,fluxes)
	call wrflx(nbbox,it,nlvdim,nsect,ivar,nlayers,fluxes)
	call fluxes_aver(1,nbc_ob,nlayers_ob,nnob,masst_ob,fluxes_ob)

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

	call boxes_3d_aver(nbox,nnn,valt)	!temp in box
	call boxes_3d_aver(nbox,nnn,vals)	!salt in box
	call boxes_2d_aver(nbox,nnn,vale)	!zeta in box
	call boxes_3d_aver(nbox,nnn,valv)	!vel speed in box

	call boxes_multi_aver(nbox,7,nnn,valm)	!meteo in box

	!call boxes_3d_aver(nbox,nnn,valw)
	call boxes_3d_aver(nbox,nnn,valvis)
	call boxes_3d_aver(nbox,nnn,valdif)

	val2d = barea*(eta_act-eta_old)
	call boxes_mass_balance_2d(dtbox,bvolume,fluxes,fluxes_ob,val2d)
	if( nlv .gt. 1  ) then
	  call boxes_mass_balance_3d(dtbox,bvolume
     +					,nblayers,nlayers
     +					,fluxes,fluxes_ob,val2d,valw)
	end if

c	-------------------------------------------------------
c	write results
c	-------------------------------------------------------

	call box_write_2d(it,fluxes,valt,vals,eta_act,valv,fluxes_ob)
	call box_write_meteo(it,valm)

	if( nlv .gt. 1 .and. bbox3d ) then
	  call box_write_3d(it,nblayers,nlayers
     +			,fluxes,valt,vals,vale,valv,fluxes_ob)
	  call box_write_vertical(it,nblayers
     +			,valw,valvis,valdif)
	end if

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

	call fluxes_init(nlvdim,nsect,nlayers,nrm,masst)
	call fluxes_init(1,nbc_ob,nlayers_ob,nnob,masst_ob)

	if( ibarcl .gt. 0 ) then
	  call fluxes_init(nlvdim,nsect,nlayers,nrs,saltt)
	  call fluxes_init(nlvdim,nsect,nlayers,nrt,tempt)
	end if

	if( iconz .eq. 1 ) then
	  call fluxes_init(nlvdim,nsect,nlayers,nrc,conzt)
	end if

	nnn = 0
	dtbox = 0.
	call boxes_3d_init(nbox,valt)	!temp in box
	call boxes_3d_init(nbox,vals)	!salt in box
	call boxes_2d_init(nbox,vale)	!zeta in box
	call boxes_3d_init(nbox,valv)	!vel speed in box

	call boxes_multi_init(nbox,7,valm)	!meteo

	call boxes_3d_init(nbox,valw)	!vertical vel
	call boxes_3d_init(nbox,valvis)	!viscosity
	call boxes_3d_init(nbox,valdif)	!diffusivity

	eta_old = eta_act

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

	use basin

	implicit none

	include 'param.h'
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

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer id,n,ib1,ib2
	integer list(n)
	integer ierr

	integer i

	if( id .gt. nscboxdim ) ierr = ierr + 1
	if( kfluxm + n + 1 .gt. nfxboxdim ) ierr = ierr + 1

	if( ierr .eq. 0 ) then	! no errors -> insert
	    isects(1,id) = n		!total number of nodes
	    isects(2,id) = kfluxm+1	!pointer into kflux
	    isects(3,id) = ib1		!from box
	    isects(4,id) = ib2		!to box
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

	use basin

	implicit none

	include 'param.h'
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

	subroutine box_make_stats(nblayers,barea,bvolume,bdepth)

c computes max layers for each box

	use mod_depth
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer nblayers(nbxdim)		!number of layers in boxes
	real barea(nbxdim)			!area of boxes
	real bvolume(nbxdim)			!volume of boxes
	real bdepth(nbxdim)			!depth of boxes


	integer ib,lmax,ie
	real area,hdep

	do ib=1,nbxdim
	  nblayers(ib) = 0
	  barea(ib) = 0.
	  bvolume(ib) = 0.
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  hdep = hev(ie)
	  area = 12.*ev(10,ie)
	  barea(ib) = barea(ib) + area
	  bvolume(ib) = bvolume(ib) + area * hdep
	  lmax = ilhv(ie)
	  nblayers(ib) = max(nblayers(ib),lmax)
	end do

	do ib=1,nbxdim
	  bdepth(ib) = bvolume(ib) / barea(ib)
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_write_stats(itanf,itend,date,time
     +					,nlayers,nblayers
     +					,barea,bvolume,bdepth)

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer itanf,itend
	integer date,time
	integer nlayers(nscboxdim)		!number of layers in section
	integer nblayers(nbxdim)		!number of layers in boxes
	real barea(nbxdim)			!area of boxes
	real bvolume(nbxdim)			!volume of boxes
	real bdepth(nbxdim)			!depth of boxes

	integer ib
	integer is,ib1,ib2
	character*20 line

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_stats.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) date,time
	call dtsgf(itanf,line)
	write(iu,*) itanf,'  ',line
	call dtsgf(itend,line)
	write(iu,*) itend,'  ',line

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,nblayers(ib),barea(ib),bvolume(ib),bdepth(ib)
	end do

        write(iu,*) nsect
        do is=1,nsect
          ib1 = isects(3,is)
          ib2 = isects(4,is)
          write(iu,*) ib1,ib2,nlayers(is),nblayers(ib1),nblayers(ib2)
        end do

	close(iu)

	end

c******************************************************************

	subroutine box_write_init(eta_act)

c writes initial conditions - still to be done : T/S

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real eta_act(nbxdim)			!water level of first time step

	integer ib
	integer is,ib1,ib2
	character*20 line

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_init.txt'
	  iu = ifileo(135,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,eta_act(ib)
	end do

	close(iu)

	end

c******************************************************************

	subroutine box_write_2d(it,fluxes,valt,vals,vale,valv,fluxes_ob)

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer it
	real fluxes(0:nlvdim,3,nscboxdim)
	real valt(0:nlvdim,nbxdim)			!temperature
	real vals(0:nlvdim,nbxdim)			!salinity
	real vale(nbxdim)				!zeta
	real valv(0:nlvdim,nbxdim)			!velocity speed
	real fluxes_ob(0:1,3,nbcdim)

	integer ib,is,ib1,ib2,ii,itype

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_2d.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) it

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,*) ib,valt(0,ib),vals(0,ib),vale(ib),valv(0,ib)
	end do

	write(iu,*) nsect
	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  write(iu,*) ib1,ib2,(fluxes(0,ii,is),ii=1,3)
	end do

	write(iu,*) nbc_ob
	do is=1,nbc_ob
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  write(iu,*) ib1,ib2,(fluxes_ob(0,ii,is),ii=1,3)
	end do

	end

c******************************************************************

	subroutine box_write_3d(it,nblayers,nlayers
     +			,fluxes,valt,vals,vale,valv,fluxes_ob)

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer it
	integer nblayers(nbxdim)		!number of layers in box
	integer nlayers(nscboxdim)		!number of layers in section
	real fluxes(0:nlvdim,3,nscboxdim)
	real valt(0:nlvdim,nbxdim)			!temperature
	real vals(0:nlvdim,nbxdim)			!salinity
	real vale(nbxdim)				!zeta
	real valv(0:nlvdim,nbxdim)			!velocity speed
	real fluxes_ob(0:1,3,nbcdim)

	integer ib,is,ib1,ib2,ii,itype
	integer lmax,l

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_3d.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) it

	write(iu,*) nbox
	do ib=1,nbox
	  lmax = nblayers(ib)
	  write(iu,*) ib,lmax,vale(ib)
	  do l=1,lmax
	    write(iu,*) l,valt(l,ib),vals(l,ib),valv(l,ib)
	  end do
	end do

	write(iu,*) nsect
	do is=1,nsect
	  lmax = nlayers(is)
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  write(iu,*) ib1,ib2,lmax
	  do l=1,lmax
	    write(iu,*) l,(fluxes(l,ii,is),ii=1,3)
	  end do
	end do

	write(iu,*) nbc_ob
	do is=1,nbc_ob
	  itype = iscbnd(2,is)
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  if( itype .eq. 1 ) then	!z-boundary, ignore
	    ib1 = 0
	    ib2 = 0
	  end if
	  write(iu,*) ib1,ib2,1
	  write(iu,*) 0,(fluxes_ob(0,ii,is),ii=1,3)
	end do

	end

c******************************************************************

	subroutine box_write_vertical(it,nblayers
     +			,valw,valvis,valdif)

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer it
	integer nblayers(nbxdim)		!number of layers in box
	real valw(0:nlvdim,nbxdim)			!temperature
	real valvis(0:nlvdim,nbxdim)			!salinity
	real valdif(0:nlvdim,nbxdim)			!velocity speed

	integer ib,is,ib1,ib2,ii,itype
	integer lmax,l,max

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_vertical.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) it

	write(iu,*) nbox
	do ib=1,nbox
	  lmax = nblayers(ib)
	  write(iu,*) ib,lmax
	  do l=1,lmax
	    write(iu,*) l,valw(l,ib),valvis(l,ib),valdif(l,ib)
	  end do
	end do

	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine box_eta(zev,val)

c computes average zeta values for box

	use mod_geom_dynamic
	use evgeom
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

        real zev(3,nel)
	real val(nbxdim)


	double precision vald(nbxdim)
	double precision vold(nbxdim)

	integer ib,ie,ii
	real vol,area3,z

	do ib=1,nbox
	  vald(ib) = 0.
	  vold(ib) = 0.
	end do

	do ie=1,nel
	  !if( iwegv(ie) .eq. 0 ) then	!only for wet elements
	    ib = iboxes(ie)
	    area3 = 4.*ev(10,ie)
	    do ii=1,3
	      z = zev(ii,ie)
	      vol = area3
	      vald(ib) = vald(ib) + vol*z
	      vold(ib) = vold(ib) + vol
	    end do
	  !end if
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

	subroutine box_vel(val)

c computes average velocity values (speed) for box

	use mod_layer_thickness
	use mod_hydro
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real val(0:nlvdim,nbxdim)


	real h,u,v
	real velspeed(nlvdim)

	double precision vald(0:nlvdim,nbxdim)
	double precision vold(0:nlvdim,nbxdim)

	integer ib,ie,ii
	integer lmax,l,k
	real vol,area3

	velspeed = 0.

	do ib=1,nbox
	  do l=0,nlvdim
	    vald(l,ib) = 0.
	    vold(l,ib) = 0.
	  end do
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
	    h = hdenv(l,ie)
	    u = utlnv(l,ie)/h
	    v = vtlnv(l,ie)/h
	    velspeed(l) = sqrt(u*u+v*v)
	  end do
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vald(0,ib) = vald(0,ib) + vol*velspeed(l)
	      vold(0,ib) = vold(0,ib) + vol
	      vald(l,ib) = vald(l,ib) + vol*velspeed(l)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	do ib=1,nbox
	  do l=0,nlvdim
	    if( vold(l,ib) .gt. 0. ) then
	      val(l,ib) = vald(l,ib)/vold(l,ib)
	    else
	      val(l,ib) = 0.
	    end if
	  end do
	end do

	end

c******************************************************************

	subroutine box_3d_vertical(scalar,val,iaver)

c computes average scalar values for box

	use mod_depth
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real scalar(0:nlvdim,nkndim)
	real val(0:nlvdim,nbxdim)
	integer iaver			!average or only accumulate


	double precision vald(0:nlvdim,nbxdim)
	double precision vold(0:nlvdim,nbxdim)

	integer ib,ie,lmax,l,ii,k
	real vol,area3

	do ib=1,nbox
	  do l=0,nlvdim
	    vald(l,ib) = 0.
	    vold(l,ib) = 0.
	  end do
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	if( iaver .le. 0 ) return	!accumulate

	do ib=1,nbox
	  do l=0,nlvdim
	    if( vold(l,ib) .gt. 0. ) then
	      val(l,ib) = vald(l,ib)/vold(l,ib)
	    else
	      val(l,ib) = 0.
	    end if
	  end do
	end do

	end

c******************************************************************

	subroutine box_3d_scalar(scalar,val)

c computes average scalar values for box

	use mod_layer_thickness
	use evgeom
	use levels
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real scalar(nlvdim,nkndim)
	real val(0:nlvdim,nbxdim)


	double precision vald(0:nlvdim,nbxdim)
	double precision vold(0:nlvdim,nbxdim)

	integer ib,ie,lmax,l,ii,k
	real vol,area3

	do ib=1,nbox
	  do l=0,nlvdim
	    vald(l,ib) = 0.
	    vold(l,ib) = 0.
	  end do
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    do l=1,lmax
	      vol = area3 * hdknv(l,k)
	      vald(0,ib) = vald(0,ib) + vol*scalar(l,k)
	      vold(0,ib) = vold(0,ib) + vol
	      vald(l,ib) = vald(l,ib) + vol*scalar(l,k)
	      vold(l,ib) = vold(l,ib) + vol
	    end do
	  end do
	end do

	do ib=1,nbox
	  do l=0,nlvdim
	    if( vold(l,ib) .gt. 0. ) then
	      val(l,ib) = vald(l,ib)/vold(l,ib)
	    else
	      val(l,ib) = 0.
	    end if
	  end do
	end do

	end

c******************************************************************

	subroutine box_2d_scalar(scalar,val)

c computes average scalar values for box

	use evgeom
	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real scalar(*)
	real val(*)

	double precision vald(nbxdim)
	double precision vold(nbxdim)

	integer ib,ie,ii,k
	real vol,area3

	do ib=1,nbox
	  vald(ib) = 0.
	  vold(ib) = 0.
	end do

	do ie=1,nel
	  ib = iboxes(ie)
	  area3 = 4.*ev(10,ie)
	  vol = area3
	  do ii=1,3
	    k = nen3v(ii,ie)
	    vald(ib) = vald(ib) + vol*scalar(k)
	    vold(ib) = vold(ib) + vol
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

	subroutine boxes_2d_init(n,val)

	implicit none

	integer n
	real val(n)

	integer i

	do i=1,n
	  val(i) = 0
	end do

	end

c******************************************************************

	subroutine boxes_2d_accum(n,val,value)

	implicit none

	integer n
	real val(n)
	real value(n)

	integer i

	do i=1,n
	  val(i) = val(i) + value(i)
	end do

	end

c******************************************************************

	subroutine boxes_2d_aver(n,na,val)

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

	subroutine boxes_3d_init(n,val)

	implicit none

	include 'param.h'

	integer n
	real val(0:nlvdim,n)

	integer i,l

	do i=1,n
	  do l=0,nlvdim
	    val(l,i) = 0
	  end do
	end do

	end

c******************************************************************

	subroutine boxes_3d_accum(n,val,value)

	implicit none

	include 'param.h'

	integer n
	real val(0:nlvdim,n)
	real value(0:nlvdim,n)

	integer i,l

	do i=1,n
	  do l=0,nlvdim
	    val(l,i) = val(l,i) + value(l,i)
	  end do
	end do

	end

c******************************************************************

	subroutine boxes_3d_aver(n,na,val)

	implicit none

	include 'param.h'

	integer n,na
	real val(0:nlvdim,n)

	integer i,l

	if( na .le. 0 ) return

	do i=1,n
	  do l=0,nlvdim
	    val(l,i) = val(l,i) / na
	  end do
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
c
c these fluxes are barotropic, so we insert then in 0 (baro) and 1 (1 layer)

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
c******************************************************************
c******************************************************************
c meteo routines
c******************************************************************
c******************************************************************
c******************************************************************

        subroutine boxes_multi_init(n,nj,val)

c as boxes_init, but for multiple values (only 2d)

        implicit none

        integer n,nj
        real val(nj,n)

        integer i,j

        do i=1,n
	  do j=1,nj
            val(j,i) = 0
	  end do
        end do

        end

c******************************************************************

	subroutine boxes_multi_aver(n,nj,na,val)

c as boxes_aver, but for multiple values (only 2d)

	implicit none

	integer n,nj,na
	real val(nj,n)

	integer i,j

	if( na .le. 0 ) return

	do i=1,n
	  do j=1,nj
            val(j,i) = val(j,i) / na
	  end do
	end do

	end

c******************************************************************

	subroutine boxes_meteo_accum(nbox,valm,val,ievap)

	use mod_meteo

	implicit none

	integer nbox
	real valm(7,nbox)
	real val(nbox)
	integer ievap

        real zconv
        parameter( zconv = 1.e-3 / 86400. )

	include 'param.h'

	integer i,ip
	real econv

	econv = 1.
	if( ievap .le. 0 ) econv = 0.

	ip = 1
	call box_2d_scalar(metrad,val)
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)
	end do

	ip = 2
	call box_2d_scalar(mettair,val)
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)
	end do

	ip = 3
	call box_2d_scalar(methum,val)
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)
	end do

	ip = 4
	call box_2d_scalar(metws,val)
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)
	end do

	ip = 5
	call box_2d_scalar(metcc,val)
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)
	end do

	ip = 6
	call box_2d_scalar(metrain,val)			![mm/day]
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)*zconv	![m/s]
	end do

	ip = 7
	call box_2d_scalar(evapv,val)			![m/s]
	do i=1,nbox
	  valm(ip,i) = valm(ip,i) + val(i)*econv
	end do

	end

c******************************************************************

	subroutine box_write_meteo(it,valm)

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	integer it
	real valm(7,nbxdim)			!meteo variables

	integer ib,j

	integer ifileo
	character*80 file
	integer iu
	save iu
	data iu /0/

	if( iu == 0 ) then
	  file = 'boxes_meteo.txt'
	  iu = ifileo(0,file,'formatted','new')
	  if( iu <= 0 ) stop 'error stop boxes: opening file'
	end if

	write(iu,*) it

	write(iu,*) nbox
	do ib=1,nbox
	  write(iu,1001) ib,(valm(j,ib),j=1,7)
	end do

 1000	format(i14,7g12.4)
 1001	format(i4,1x,5f8.3,2g12.4)
	end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine boxes_mass_balance_2d(dtbox,bvolume
     +					,fluxes,fluxes_ob,vale)

c computes mass balance

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real dtbox
	real bvolume(nbxdim)
	real fluxes(0:nlvdim,3,nscboxdim)
	real fluxes_ob(0:1,3,nbcdim)
	real vale(nbxdim)		!volume difference due to eta

	integer ib,is,ib1,ib2
	real volf,vole,volb,vola,vold,volr
	real errmax
	double precision vol(nbxdim)
	double precision flux
	logical bdebug

	real eps
	parameter (eps=1.E-5)

	bdebug = .true.
	bdebug = .false.

	do ib=1,nbox
	  vol(ib) = 0.
	end do

	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  flux = dtbox * fluxes(0,1,is)
	  if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxes_ob(0,1,is)
	  if( ib1 > 0 ) vol(ib1) = vol(ib1) - flux
	  if( ib2 > 0 ) vol(ib2) = vol(ib2) + flux
	end do

	errmax = 0.
	do ib=1,nbox
	  volf = vol(ib)
	  vole = vale(ib)
	  volb = bvolume(ib)
	  vola = max(abs(volf),abs(vole))
	  vold = abs(volf-vole)
	  volr = 0.
	  !volr = vold/vola
	  volr = vold/volb	!use total volume
	  if( volr .gt. eps ) then
	    write(6,*) '*** mass balance error in box ',ib
	    write(6,*) '***',ib,volf,vole,vold,volr
	  end if
	  errmax = max(errmax,abs(volr))
	  !write(115,*) ib,volf,vole,vold,volr
	end do

	if( bdebug ) write(6,*) 'errmax 2d: ',errmax
	if( errmax .gt. eps ) then
	  write(6,*) 'errmax 2d: ',errmax
	end if

	end

c******************************************************************

	subroutine boxes_mass_balance_3d(dtbox,bvolume
     +					,nblayers,nlayers
     +					,fluxes,fluxes_ob
     +					,vale,wflux)

c computes mass balance

	use basin

	implicit none

	include 'param.h'
	include 'subboxa.h'

	real dtbox
	real bvolume(nbxdim)
	integer nblayers(nbxdim)
	integer nlayers(nscboxdim)		!number of layers in section
	real fluxes(0:nlvdim,3,nscboxdim)
	real fluxes_ob(0:1,3,nbcdim)
	real vale(nbxdim)			!volume difference due to eta
	real wflux(0:nlvdim,nbxdim)		!vertical fluxes [m**3/s]

	integer ib,is,ib1,ib2
	integer l,lmax,ltot
	real volf,vole,volb,vola,vold,volr,volv
	real wtop,errmax
	double precision vol(nlvdim,nbxdim)
	double precision flux

	logical bdebug

	real eps
	parameter (eps=1.E-5)

	bdebug = .true.
	bdebug = .false.

	ltot = 0
	do ib=1,nbox
	  ltot = max(ltot,nblayers(ib))
	end do

	do ib=1,nbox
	  do l=1,ltot
	    vol(l,ib) = 0.
	  end do
	end do

	do is=1,nsect
	  ib1 = isects(3,is)
	  ib2 = isects(4,is)
	  lmax = nlayers(is)
	  do l=1,lmax
	    flux = dtbox * fluxes(l,1,is)
	    if( ib1 > 0 ) vol(l,ib1) = vol(l,ib1) - flux
	    if( ib2 > 0 ) vol(l,ib2) = vol(l,ib2) + flux
	    if( bdebug .and. abs(flux) > 0. ) then
	      if( nblayers(ib1) < l .or. nblayers(ib2) < l ) then
		write(6,'(a,6i5,g14.4)') 'non existing flux: '
     +				,is,ib1,ib2
     +				,nblayers(ib1),nblayers(ib2),l,flux
	      end if
	    end if
	  end do
	end do

	do is=1,nbc_ob
	  ib1 = iscbnd(3,is)
	  ib2 = iscbnd(4,is)
	  flux = dtbox * fluxes_ob(1,1,is)
	  if( ib1 > 0 ) vol(1,ib1) = vol(1,ib1) - flux
	  if( ib2 > 0 ) vol(1,ib2) = vol(1,ib2) + flux
	end do

	do ib=1,nbox
	  wtop = 0.
	  wflux(ltot,ib) = wtop		!bottom
	  do l=ltot,1,-1
	    wtop = wtop + vol(l,ib)
	    wflux(l-1,ib) = -wtop / dtbox
	  end do
	end do

	errmax = 0.
	do ib=1,nbox
	  lmax = nblayers(ib)
	  do l=lmax+1,ltot
	    if( bdebug .and. vol(l,ib) /= 0. ) then
	      write(6,'(a,3i5,g14.4)') 'non existing box: '
     +				,ib,lmax,l,vol(l,ib)
	    end if
	  end do
	  
	  volb = bvolume(ib)
	  do l=1,ltot
	    volf = vol(l,ib)
	    volv = dtbox * (wflux(l-1,ib) - wflux(l,ib))
	    vold = abs(volf + volv)
	    volr = vold/volb
	    !write(116,'(2i5,4g14.4)') ib,l,volf,volv,vold,volr
	    errmax = max(errmax,abs(volr))
	  end do

	  wtop = -dtbox * wflux(0,ib)
	  vole = vale(ib)
	  vold = abs(wtop - vole)
	  volr = vold/volb
	  !write(116,'(a,i5,4g14.4)') 'first layer: ',ib,wtop,vole,vold,volr
	  errmax = max(errmax,abs(volr))

	end do

	if( bdebug ) write(6,*) 'errmax 3d: ',errmax
	if( errmax .gt. eps ) then
	  write(6,*) 'errmax 3d: ',errmax
	end if

	end

c******************************************************************

