c
c $Id: subcon1.f,v 1.21 2010-03-11 15:36:39 georg Exp $
c
c routines for concentration (utilities) (old newcon1.f)
c
c contents :
c
c subroutine conini0(nlvdi,c,cref)                      sets initial conditions
c subroutine conini(nlvdi,c,cref,cstrat)		sets initial conditions
c
c subroutine conbnd(nlvdi,c,rbc)			boundary condition
c subroutine con3bnd(nlvdi,c,nlvbnd,rbc)		boundary condition (3D)
c
c subroutine confop(iu,itmcon,idtcon,nlv,nvar,type)	opens (NOS) file
c subroutine confil(iu,itmcon,idtcon,ivar,nlvdi,c)	writes NOS file
c subroutine conwrite(iu,type,nvar,ivar,nlvdi,c)        shell for writing file
c
c subroutine conmima(nlvdi,c,cmin,cmax)                 computes min/max
c subroutine conmimas(nlvdi,c,cmin,cmax)                computes scalar min/max
c
c notes :
c
c conbnd and con3bnd are not used and can be deleted
c
c revision log :
c
c 19.08.1998	ggu	call to conzfi changed
c 20.08.1998	ggu	makew removed (routine used is sp256w)
c 24.08.1998    ggu     levdbg used for debug
c 26.08.1998    ggu     subroutine convol, tstvol transferred to newchk
c 26.08.1998    ggu     all subroutines re-written more generally
c 26.01.1999    ggu     can be used also with 2D routines
c 16.11.2001    ggu     subroutine conmima and diffstab
c 05.12.2001    ggu     new routines diffstab,diffstab1,difflimit
c 11.10.2002    ggu     commented diffset
c 09.09.2003    ggu     new routine con3bnd
c 10.03.2004    ggu     new routine conwrite()
c 13.03.2004    ggu     new routines set_c_bound, distribute_vertically
c 13.03.2004    ggu     exec routine con3bnd() only for level BC (LEVELBC)
c 14.03.2004    ggu     new routines open_b_flux
c 05.01.2005    ggu     routine to write 2d nos file into subnosa.f
c 07.01.2005    ggu     routine diffwrite deleted
c 14.01.2005    ggu     new file for diffusion routines (copied to subdif.f)
c 23.03.2006    ggu     changed time step to real
c 31.05.2007    ggu     reset BC of flux type to old way (DEBHELP)
c 07.04.2008    ggu     deleted set_c_bound
c 08.04.2008    ggu     cleaned, deleted distribute_vertically, open_b_flux
c 09.10.2008    ggu&ccf call to confop changed -> nlv
c 20.11.2009    ggu	in conwrite only write needed (nlv) layers
c
c*****************************************************************

	subroutine conini0(nlvdi,c,cref)

c sets initial conditions (no stratification)

	implicit none

	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!variable to initialize
	real cref		!reference value
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
c local
	integer k,l
	real depth,hlayer

	do k=1,nkn
	  do l=1,nlvdi
	    c(l,k) = cref
	  end do
	end do

	end

c*****************************************************************

	subroutine conini(nlvdi,c,cref,cstrat,hdko)

c sets initial conditions (with stratification)

	implicit none

	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!variable to initialize
	real cref		!reference value
	real cstrat		!stratification [conc/km]
	real hdko(nlvdi,1)	!layer thickness
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
c local
	integer k,l
	real depth,hlayer

	do k=1,nkn
	  depth=0.
	  do l=1,nlvdi
	    hlayer = 0.5 * hdko(l,k)
	    depth = depth + hlayer
	    c(l,k) = cref + cstrat*depth/1000.
	    depth = depth + hlayer
	  end do
	end do

	end

c*************************************************************

	subroutine conbnd(nlvdi,c,rbc)

c implements boundary condition (simplicistic version)

	implicit none

c arguments
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!concentration (cconz,salt,temp,...)
	real rbc(1)		!boundary condition
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
        integer ilhkv(1)
        common /ilhkv/ilhkv
c local
	integer k,l,lmax
	real rb

	do k=1,nkn
	  if( rbc(k) .ne. flag ) then
	    rb = rbc(k)
	    lmax=ilhkv(k)
	    !write(6,*) 'conbnd: ',k,lmax,rb
	    do l=1,lmax
		c(l,k) = rb
	    end do
	  end if
	end do

	end

c*************************************************************

	subroutine con3bnd(nlvdi,c,nlvbnd,rbc)

c implements boundary condition (simplicistic 3D version)

	implicit none

c arguments
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!concentration (cconz,salt,temp,...)
	integer nlvbnd		!vertical dimension of boundary conditions
	real rbc(nlvbnd,1)	!boundary condition
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real rzv(1)
        common /rzv/rzv
c local
	integer k,l,lmax
	real rb
        integer ipext

	if( nlvbnd .ne. 1 .and. nlvbnd .ne. nlvdi ) then
	  write(6,*) 'nlvdi,nlvbnd: ',nlvdi,nlvbnd
	  stop 'error stop con3bnd: impossible nlvbnd'
	end if
	if( nlvbnd .ne. nlvdi ) then
	  write(6,*) 'nlvdi,nlvbnd: ',nlvdi,nlvbnd
	  stop 'error stop con3bnd: only 3D boundary conditions'
	end if

	do k=1,nkn
	 if( rzv(k) .ne. flag ) then    !only level BC  !LEVELBC !DEBHELP
	  if( rbc(1,k) .ne. flag ) then
	    lmax=ilhkv(k)
            !write(94,*) 'con3bnd: ',k,ipext(k),lmax,nlvbnd,rbc(1,k)
	    if( nlvbnd .eq. 1 ) then
	      rb = rbc(1,k)
	      do l=1,lmax
		c(l,k) = rb
	      end do
	    else
	      do l=1,lmax
		c(l,k) = rbc(l,k)
	      end do
	    end if
	  end if
	 end if
	end do

	end

c*************************************************************

	subroutine confop(iu,itmcon,idtcon,nlv,nvar,type)

c opens (NOS) file

c on return iu = -1 means that no file has been opened and is not written

	implicit none

	integer iu		!unit				       (in/out)
	integer itmcon		!time of first write		       (in/out)
	integer idtcon		!time intervall of writes	       (in/out)
	integer nlv		!vertical dimension of scalar          (in)
	integer nvar		!total number of variables to write    (in)
	character*(*) type	!extension of file		       (in)

        character*80 descrp
        common /descrp/ descrp
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer ilhkv(1)
        common /ilhkv/ilhkv
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev

	integer ierr
	character*80 dir,nam,file

	integer ifileo

c check idtcon and itmcon and adjust

        if( idtcon .le. 0 ) iu = -1
        if( itmcon+idtcon .gt. itend ) iu = -1

        if( iu .eq. -1 ) return

        if( itmcon .lt. itanf ) itmcon = itanf
	if( idtcon .lt. idt ) idtcon = idt

c construct file name

c        call getfnm('datdir',dir)
c        call getfnm('runnam',nam)
c        call mkname(dir,nam,type,file)

	call deffile(type,file)

c open file

	iu = ifileo(iu,file,'unformatted','new')
	if( iu .le. 0 ) goto 98

c write header of file

	call wfnos(iu,3,nkn,nel,nlv,nvar,descrp,ierr)
        if(ierr.gt.0) goto 99
	call wsnos(iu,ilhkv,hlv,hev,ierr)
        if(ierr.gt.0) goto 99

c write informational message to terminal

        write(6,*) 'confop: ',type,' file opened ',it

c end

	return
   98	continue
	write(6,*) 'error opening file ',file
	stop 'error stop confop'
   99	continue
	write(6,*) 'error ',ierr,' writing file ',file
	stop 'error stop confop'
	end

c*************************************************************

	subroutine confil(iu,itmcon,idtcon,ivar,nlvdi,c)

c writes NOS file

	implicit none

	integer iu		!unit
	integer itmcon		!time of first write
	integer idtcon		!time intervall of writes
	integer ivar		!id of variable to be written
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!concentration to write

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer ierr

c check if files has to be written

	if( iu .le. 0 ) return
	if( it .lt. itmcon ) return
	if( mod(it-itmcon,idtcon) .ne. 0 ) return

c write file

	call wrnos(iu,it,ivar,nlvdi,ilhkv,c,ierr)
	if(ierr.gt.0) goto 99

c write informational message to terminal

        write(6,*) 'confil: variable ',ivar,' written at ',it

c end

	return
   99	continue
	write(6,*) 'error ',ierr,' writing file at unit ',iu
	stop 'error stop confil'
	end

c*************************************************************

	subroutine conwrite(iu,type,nvar,ivar,nlvdim,c)

c shell for writing file unconditionally to disk

        implicit none

	integer iu		!unit (0 for first call, set on return)
        character*(*) type      !type of file
	integer nvar		!total number of variables
	integer ivar		!id of variable to be written
	integer nlvdim		!vertical dimension of c
	real c(nlvdim,1)	!concentration to write

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        integer itmcon,idtcon,lmax

        !if(nlvdim.ne.nlvdi) stop 'error stop conwrite: level dimension'

        itmcon = itanf
	idtcon = idt
	idtcon = 1
	lmax = min(nlvdim,nlv)

        if( iu .eq. 0 ) then
	  call confop(iu,itmcon,idtcon,lmax,nvar,type)
        end if

	call confil(iu,itmcon,idtcon,ivar,nlvdim,c)

        end

c*************************************************************

        subroutine conmima(nlvdi,c,cmin,cmax)

c computes min/max for scalar field

	implicit none

c arguments
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv
c local
	integer k,l,lmax
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
	  end do
	end do

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

        subroutine conmimas(nlvdi,c,cmin,cmax)

c computes min/max for scalar field -> writes some info

	implicit none

c arguments
	integer nlvdi		!vertical dimension of c
	real c(nlvdi,1)		!concentration (cconz,salt,temp,...)
        real cmin,cmax
c common
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it
        integer ilhkv(1)
        common /ilhkv/ilhkv
c local
	integer k,l,lmax
        integer ntot
	real cc
        logical debug
        integer kcmin,lcmin,kcmax,lcmax

        debug = .false.
        cmin = c(1,1)
        cmax = c(1,1)

        ntot = 0
	do k=1,nkn
	  lmax=ilhkv(k)
	  do l=1,lmax
	    cc = c(l,k)
            if( debug ) then
              if( cc .lt. cmin ) then
                    kcmin = k
                    lcmin = l
              end if
              if( cc .gt. cmax ) then
                    kcmax = k
                    lcmax = l
              end if
            end if
            cmin = min(cmin,cc)
            cmax = max(cmax,cc)
            if( cc .le. 0. ) then
                    ntot = ntot + 1
                    write(96,*) it,l,k,cc,ntot
            end if
	  end do
	end do

        if( ntot .gt. 0 ) then
                write(96,*) 'ntot: ',it,ntot
        end if

        if( debug ) then
          write(6,*) 'conmima: ',kcmin,lcmin,cmin
          write(6,*) 'conmima: ',kcmax,lcmax,cmax
        end if

        end

c*************************************************************

