!
! $Id: ousutil.f,v 1.3 2009-09-14 08:20:58 georg Exp $
!
! utilities for OUS files
!
! revision log :
!
! 16.12.2010	ggu	copied from ousextr_gis.f
! 03.06.2011	ggu	some routines transfered to genutil.f
! 08.06.2011	ggu	new routine transp2nodes()
! 10.11.2011    ggu     new routines for hybrid levels
! 02.12.2011    ggu     bug fix for call to get_sigma_info() (missing argument)
! 21.01.2013    ggu     added two new routines comp_vel2d, comp_barotropic
! 05.09.2013    ggu     new call to get_layer_thickness()
! 20.01.2014    ggu     new helper routines
! 29.04.2015    ggu     new helper routines for start/end time in file
! 11.09.2015    ggu     weight and hl are local, new velzeta2scal()
!
!******************************************************************
!------------------------------------------------------------------------------
        module ousutil
!------------------------------------------------------------------------------
        contains
!------------------------------------------------------------------------------

        subroutine transp2vel(nel,nkn,nlv,nlvddi,hev,zenv,nen3v,ilhv,hlv,utlnv,vtlnv,uprv,vprv)

! transforms transports at elements to velocities at nodes

        use evgeom_2nc
        use sigma
        use output_util

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvddi
        double precision hev(1)
        double precision zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	double precision hlv(1)
        double precision utlnv(nlvddi,1)
        double precision vtlnv(nlvddi,1)
        double precision uprv(nlvddi,1)
        double precision vprv(nlvddi,1)

        double precision weight(nlvddi,nkn)		!aux variable for weights
	double precision hl(nlvddi)			!aux variable for real level thickness

	logical bsigma
        integer ie,ii,k,l,lmax,nsigma,nlvaux
        double precision hmed,u,v,area,zeta
	double precision hsigma

	call get_sigma_info(nlvaux,nsigma,hsigma)
	if( nlvaux .gt. nlvddi ) stop 'error stop transp2vel: nlvddi'
	bsigma = nsigma .gt. 0

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    uprv(l,k) = 0.d0
	    vprv(l,k) = 0.d0
	  end do
	end do
	      
        do ie=1,nel

	  area = area_elem(ie)
	  lmax = ilhv(ie)
	  call compute_levels_on_element(ie,zenv,zeta)
	  call get_layer_thickness(lmax,nsigma,hsigma,zeta,hev(ie),hlv,hl)
	  !call get_layer_thickness_e(ie,lmax,bzeta,nsigma,hsigma,hl)

	  do l=1,lmax
	    hmed = hl(l)
	    u = utlnv(l,ie) / hmed
	    v = vtlnv(l,ie) / hmed
	    do ii=1,3
	      k = nen3v(ii,ie)
	      uprv(l,k) = uprv(l,k) + area * u
	      vprv(l,k) = vprv(l,k) + area * v
	      weight(l,k) = weight(l,k) + area
	    end do
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    area = weight(l,k)
	    if( area .gt. 0. ) then
	      uprv(l,k) = uprv(l,k) / area
	      vprv(l,k) = vprv(l,k) / area
	    end if
	  end do
	end do
	      
	end

!******************************************************************

        subroutine transp2nodes(nel,nkn,nlv,nlvddi,hev,zenv,nen3v,ilhv,hlv,utlnv,vtlnv,utprv,vtprv,weight)

! transforms transports at elements to transports at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvddi
        double precision hev(1)
        double precision zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	double precision hlv(1)
        double precision utlnv(nlvddi,1)
        double precision vtlnv(nlvddi,1)
        double precision utprv(nlvddi,1)
        double precision vtprv(nlvddi,1)
        double precision weight(nlvddi,1)

        integer ie,ii,k,l,lmax
        double precision u,v,w

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    utprv(l,k) = 0.
	    vtprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    u = utlnv(l,ie)
	    v = vtlnv(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      utprv(l,k) = utprv(l,k) + u
	      vtprv(l,k) = vtprv(l,k) + v
	      weight(l,k) = weight(l,k) + 1.
	    end do
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    w = weight(l,k)
	    if( w .gt. 0. ) then
	      utprv(l,k) = utprv(l,k) / w
	      vtprv(l,k) = vtprv(l,k) / w
	    end if
	  end do
	end do
	      
	end

!***************************************************************

        subroutine comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v,umin,vmin,umax,vmax)

! computes 2D velocities from 2D transports - returns result in u2v,v2v

        implicit none

        integer nel
        double precision hev(1)
        double precision zenv(3,1)
        double precision ut2v(1)
        double precision vt2v(1)
	double precision u2v(1), v2v(1)
        double precision umin,vmin
        double precision umax,vmax

        integer ie,ii
        double precision zmed,hmed,u,v

	umin = +1.e+30
	vmin = +1.e+30
	umax = -1.e+30
	vmax = -1.e+30

        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.
          hmed = hev(ie) + zmed

          u = ut2v(ie) / hmed
          v = vt2v(ie) / hmed

	  u2v(ie) = u
	  v2v(ie) = v

          umin = min(umin,u)
          vmin = min(vmin,v)
          umax = max(umax,u)
          vmax = max(vmax,v)
        end do

        end

!***************************************************************

	subroutine comp_barotropic(nel,nlvddi,ilhv,utlnv,vtlnv,ut2v,vt2v)

! computes barotropic transport

	implicit none

	integer nel,nlvddi
	integer ilhv(1)
	double precision utlnv(nlvddi,1)
	double precision vtlnv(nlvddi,1)
	double precision ut2v(1)
	double precision vt2v(1)

	integer ie,l,lmax
	double precision utot,vtot

	do ie=1,nel
	  lmax = ilhv(ie)
	  utot = 0.
	  vtot = 0.
	  do l=1,lmax
	    utot = utot + utlnv(l,ie)
	    vtot = vtot + vtlnv(l,ie)
	  end do
	  ut2v(ie) = utot
	  vt2v(ie) = vtot
	end do

	end

!***************************************************************

	subroutine compute_volume(nel,zenv,hev,volume)

	use evgeom_2nc

	implicit none

	integer nel
	double precision zenv(3,1)
	double precision hev(1)
	double precision volume

	integer ie,ii
	double precision zav,area
	double precision vol,voltot,areatot

	voltot = 0.
	areatot = 0.

	do ie=1,nel
	  zav = 0.
	  do ii=1,3
	    zav = zav + zenv(ii,ie)
	  end do
	  area = area_elem(ie)
	  vol = area * (hev(ie) + zav/3.)
	  voltot = voltot + vol
	  !areatot = areatot + area
	end do

	volume = voltot

	end

!***************************************************************

	subroutine compute_volume_ia(ia,zenv,volume,areaall)

	use evgeom_2nc
	use basin

	implicit none

	integer ia		!area code for which not to compute volume
	double precision zenv(3,nel)
	double precision volume,areaall

	integer ie,ii
	double precision h,area
	double precision vol,voltot,areatot

	voltot = 0.
	areatot = 0.

	do ie=1,nel
	  if( iarv(ie) == ia ) cycle
	  h = 0.
	  do ii=1,3
	    h = h + hm3v(ii,ie) + zenv(ii,ie)
	  end do
	  h = h / 3.
	  area = area_elem(ie)
	  vol = area * h
	  voltot = voltot + vol
	  areatot = areatot + area
	end do

	volume = voltot
	areaall = areatot

	end

!***************************************************************

        subroutine debug_write_node(ks,it,nrec,nknddi,nelddi,nlvddi,nkn,nel,nlv,nen3v,zenv,znv,utlnv,vtlnv)

! debug write

        implicit none

	integer ks	!internal node number to output (0 for none)
        integer it,nrec
        integer nknddi,nelddi,nlvddi,nkn,nel,nlv
        integer nen3v(3,nelddi)
        double precision znv(nknddi)
        double precision zenv(3,nelddi)
        double precision utlnv(nlvddi,nelddi)
        double precision vtlnv(nlvddi,nelddi)

        integer ie,ii,k,l
        logical bk

	if( ks .le. 0 ) return

        write(66,*) 'time: ',it,nrec
        write(66,*) 'kkk: ',ks,znv(ks)

        do ie=1,nel
          bk = .false.
          do ii=1,3
            k = nen3v(ii,ie)
            if( k .eq. ks ) then
              write(66,*) 'ii: ',ii,ie,zenv(ii,ie)
              bk = .true.
            end if
          end do
          if( bk ) then
          do l=1,nlv
            write(66,*) 'ie: ',ie,l,utlnv(l,ie),vtlnv(l,ie)
          end do
          end if
        end do

        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine write_ous_header(iu,ilhv,hlv,hev)

! other variables are stored internally
!
! must have been initialized with ous_init
! all other variables must have already been stored internally (title,date..)

        use ioous

        implicit none

        integer iu
        integer ilhv(1)
        double precision hlv(1)
        double precision hev(1)

        integer nkn,nel,nlv
        integer ierr

        call ous_get_params(iu,nkn,nel,nlv)
        call ous_write_header(iu,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) goto 99
        call ous_write_header2(iu,ilhv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 99

        return
   99   continue
        write(6,*) 'error in writing header of OUS file'
        stop 'error stop write_ous_header: writing header'
        end

!***************************************************************

        subroutine peek_ous_header(iu,nkn,nel,nlv)

! get size of data

        use ioous

        implicit none

        integer iu
        integer nkn,nel,nlv

        integer nvers
        integer ierr

        nvers = 2
        call ous_init(iu,nvers)

        call ous_read_header(iu,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) goto 99

        call ous_close(iu)
        rewind(iu)

        return
   99   continue
        write(6,*) 'error in reading header of OUS file'
        stop 'error stop peek_ous_header: reading header'
        end

!***************************************************************

        subroutine read_ous_header(iu,nknddi,nelddi,nlvddi,ilhv,hlv,hev)

! other variables are stored internally

        use ioous

        implicit none

        integer iu
        integer nknddi,nelddi,nlvddi
        integer ilhv(nknddi)
        double precision hlv(nlvddi)
        double precision hev(nelddi)

        integer nvers
        integer nkn,nel,nlv,nvar
        integer ierr
        integer l
        integer date,time
	double precision href,hzmin
        character*50 title,femver

        nvers = 2

        call ous_init(iu,nvers)

        call ous_read_header(iu,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) goto 99

        call dimous(iu,nknddi,nelddi,nlvddi)
	!call infoous(iu,6)

        call getous(iu,nvers,nkn,nel,nlv)
        call ous_get_date(iu,date,time)
        call ous_get_title(iu,title)
        call ous_get_femver(iu,femver)
        call ous_get_hparams(iu,href,hzmin)

        write(6,*) 'nvers      : ',nvers
        write(6,*) 'nkn,nel    : ',nkn,nel
        write(6,*) 'nlv        : ',nlv
        write(6,*) 'title      : ',title
        write(6,*) 'femver     : ',femver
        write(6,*) 'date,time  : ',date,time
        write(6,*) 'href,hzmin : ',href,hzmin

        call ous_read_header2(iu,ilhv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 99

        write(6,*) 'Available levels: ',nlv
        write(6,'(5g14.6)') (hlv(l),l=1,nlv)

        return
   99   continue
        write(6,*) 'error in reading header of OUS file'
        stop 'error stop read_ous_header: reading header'
        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine open_ous_type(type,status,nunit)

! open OUS file with default simulation name and given extension

        use defnames
        use fil

        implicit none

        character*(*) type,status
        integer nunit

        integer nb
        character*80 file

        call def_make(type,file)
        nb = ifileo(0,file,'unform',status)

        if( nb .le. 0 ) then
          write(6,*) 'file: ',file
          stop 'error stop open_ous_type: opening file'
        end if

        nunit = nb

        end

!***************************************************************

        subroutine open_ous_file(name,status,nunit)

        use defnames
        use fil

        implicit none

        character*(*) name,status
        integer nunit

        integer nb
        character*80 file

        call mkname(' ',name,'.ous',file)
        nb = ifileo(0,file,'unform',status)
        if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_ous_file: opening file'
	end if

        nunit = nb

        end

!***************************************************************

        subroutine qopen_ous_file(text,status,nunit)

! asks for name and opens ous file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_ous_file(name,status,nunit)

        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine ous_get_it_start(file,itstart)

! gets it of first record

        use ioous

        implicit none

        character*(*) file
        integer itstart

        integer nunit,nvers
        integer it,ierr
        character*80 title

        nvers = 2
        itstart = -1

        call open_ous_file(file,'old',nunit)
        call ous_init(nunit,nvers)
        call ous_skip_header(nunit,ierr)
        if( ierr .ne. 0 ) return
        call ous_skip_record(nunit,it,ierr)
        if( ierr .ne. 0 ) return
        itstart = it

        end

!***************************************************************

        subroutine ous_get_it_end(file,itend)

! gets it of last record

        use ioous

        implicit none

        character*(*) file
        integer itend

        integer nunit,nvers
        integer it,itlast,ierr
        character*80 title

        nvers = 2
        itend = -1
        itlast = -1

        call open_ous_file(file,'old',nunit)
        call ous_init(nunit,nvers)
        call ous_skip_header(nunit,ierr)
        if( ierr .ne. 0 ) return

    1   continue
        call ous_skip_record(nunit,it,ierr)
        if( ierr .gt. 0 ) return
        if( ierr .lt. 0 ) goto 2
        itlast = it
        goto 1
    2   continue
        itend = itlast

        end

!***************************************************************

        function check_ous_file(file)

        use fil
        use ioous

        implicit none

        logical check_ous_file
        character*(*) file

        integer nb,nvers

        check_ous_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
        call ous_is_ous_file(nb,nvers)
        close(nb)

        check_ous_file = nvers > 0
        
        end

!***************************************************************

        subroutine velzeta2scal(nel,nkn,nlv,nlvddi,nen3v,ilhkv,zenv,uprv,vprv,zv,sv,dv)

        use geometrical

	implicit none

	integer nkn,nel,nlv,nlvddi
	integer nen3v(3,nel)
	integer ilhkv(nkn)
	double precision zenv(3,nel)
	double precision uprv(nlvddi,nkn)
	double precision vprv(nlvddi,nkn)
	double precision zv(nkn)
	double precision sv(nlvddi,nkn)
	double precision dv(nlvddi,nkn)

	integer ie,k,ii,l,lmax
	double precision u,v,s,d

	zv = 1.e+30
	sv = 0.
	dv = 0.

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zv(k) = min(zv(k),zenv(ii,ie))
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    u = uprv(l,k)
	    v = vprv(l,k)
	    call c2p(u,v,s,d)	!d is meteo convention
	    d = d + 180.
	    if( d > 360. ) d = d - 360.
	    sv(l,k) = s
	    dv(l,k) = d
	  end do
	end do
	
	end

!***************************************************************

!------------------------------------------------------------------------------
        end module ousutil
!------------------------------------------------------------------------------
