
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2017,2019  Georg Umgiesser
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

c utilities for OUS files
c
c revision log :
c
c 16.12.2010	ggu	copied from ousextr_gis.f
c 03.06.2011	ggu	some routines transfered to genutil.f
c 08.06.2011	ggu	new routine transp2nodes()
c 10.11.2011	ggu	new routines for hybrid levels
c 02.12.2011	ggu	bug fix for call to get_sigma_info() (missing argument)
c 09.12.2011	ggu	changed VERS_6_1_38
c 27.01.2012	ggu	changed VERS_6_1_43
c 21.01.2013	ggu	added two new routines comp_vel2d, comp_barotropic
c 13.06.2013	ggu	changed VERS_6_1_65
c 05.09.2013	ggu	new call to get_layer_thickness()
c 12.09.2013	ggu	changed VERS_6_1_67
c 20.01.2014	ggu	new helper routines
c 28.01.2014	ggu	changed VERS_6_1_71
c 29.04.2015	ggu	new helper routines for start/end time in file
c 05.06.2015	ggu	changed VERS_7_1_12
c 10.07.2015	ggu	changed VERS_7_1_50
c 17.07.2015	ggu	changed VERS_7_1_53
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 30.07.2015	ggu	changed VERS_7_1_83
c 11.09.2015	ggu	weight and hl are local, new velzeta2scal()
c 05.11.2015	ggu	changed VERS_7_3_12
c 19.02.2016	ggu	changed VERS_7_5_2
c 19.02.2016	ggu	changed VERS_7_5_3
c 09.09.2016	ggu	changed VERS_7_5_17
c 14.11.2017	ggu	changed VERS_7_5_36
c 16.02.2019	ggu	changed VERS_7_5_60
c
c******************************************************************

        subroutine transp2vel(nel,nkn,nlv,nlvddi,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv)

c transforms transports at elements to velocities at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvddi
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvddi,1)
        real vtlnv(nlvddi,1)
        real uprv(nlvddi,1)
        real vprv(nlvddi,1)

        real weight(nlvddi,nkn)		!aux variable for weights
	real hl(nlvddi)			!aux variable for real level thickness

	logical bsigma
        integer ie,ii,k,l,lmax,nsigma,nlvaux
        real hmed,u,v,area,zeta
	real hsigma

	real area_elem

	call get_sigma_info(nlvaux,nsigma,hsigma)
	if( nlvaux .gt. nlvddi ) stop 'error stop transp2vel: nlvddi'
	bsigma = nsigma .gt. 0

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    uprv(l,k) = 0.
	    vprv(l,k) = 0.
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

c******************************************************************

        subroutine transp2nodes(nel,nkn,nlv,nlvddi,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,utprv,vtprv,weight)

c transforms transports at elements to transports at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvddi
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvddi,1)
        real vtlnv(nlvddi,1)
        real utprv(nlvddi,1)
        real vtprv(nlvddi,1)
        real weight(nlvddi,1)

        integer ie,ii,k,l,lmax
        real u,v,w

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

c***************************************************************

        subroutine comp_vel2d(nel,hev,zenv,ut2v,vt2v,u2v,v2v
     +				,umin,vmin,umax,vmax)

c computes 2D velocities from 2D transports - returns result in u2v,v2v

        implicit none

        integer nel
        real hev(1)
        real zenv(3,1)
        real ut2v(1)
        real vt2v(1)
	real u2v(1), v2v(1)
        real umin,vmin
        real umax,vmax

        integer ie,ii
        real zmed,hmed,u,v

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

c***************************************************************

	subroutine comp_barotropic(nel,nlvddi,ilhv
     +			,utlnv,vtlnv,ut2v,vt2v)

c computes barotropic transport

	implicit none

	integer nel,nlvddi
	integer ilhv(1)
	real utlnv(nlvddi,1)
	real vtlnv(nlvddi,1)
	real ut2v(1)
	real vt2v(1)

	integer ie,l,lmax
	real utot,vtot

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

c***************************************************************

	subroutine compute_volume(nel,zenv,hev,volume)

	use evgeom

	implicit none

	integer nel
	real zenv(3,1)
	real hev(1)
	real volume

	integer ie,ii
	real zav,area
	double precision vol,voltot,areatot

	real area_elem

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

c***************************************************************

	subroutine compute_volume_ia(ia,zenv,volume,areaall)

	use evgeom
	use basin

	implicit none

	integer ia		!area code for which not to compute volume
	real zenv(3,nel)
	real volume,areaall

	integer ie,ii
	real h,area
	double precision vol,voltot,areatot

	real area_elem

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

c***************************************************************

        subroutine debug_write_node(ks,it,nrec
     +		,nknddi,nelddi,nlvddi,nkn,nel,nlv
     +          ,nen3v,zenv,znv,utlnv,vtlnv)

c debug write

        implicit none

	integer ks	!internal node number to output (0 for none)
        integer it,nrec
        integer nknddi,nelddi,nlvddi,nkn,nel,nlv
        integer nen3v(3,nelddi)
        real znv(nknddi)
        real zenv(3,nelddi)
        real utlnv(nlvddi,nelddi)
        real vtlnv(nlvddi,nelddi)

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

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine write_ous_header(iu,ilhv,hlv,hev)

c other variables are stored internally
c
c must have been initialized with ous_init
c all other variables must have already been stored internally (title,date..)

        implicit none

        integer iu
        integer ilhv(1)
        real hlv(1)
        real hev(1)

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

c***************************************************************

        subroutine peek_ous_header(iu,nkn,nel,nlv)

c get size of data

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

c***************************************************************

        subroutine read_ous_header(iu,nknddi,nelddi,nlvddi,ilhv,hlv,hev)

c other variables are stored internally

        implicit none

        integer iu
        integer nknddi,nelddi,nlvddi
        integer ilhv(nknddi)
        real hlv(nlvddi)
        real hev(nelddi)

        integer nvers
        integer nkn,nel,nlv,nvar
        integer ierr
        integer l
        integer date,time
	real href,hzmin
        character*80 title,femver

        nvers = 2

        call ous_init(iu,nvers)

        call ous_read_header(iu,nkn,nel,nlv,ierr)
        if( ierr .ne. 0 ) goto 99

        call dimous(iu,nknddi,nelddi,nlvddi)

        call ous_read_header2(iu,ilhv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 99

        !write(6,*) 'levels: '
        !write(6,'(5g14.6)') (hlv(l),l=1,nlv)

        return
   99   continue
        write(6,*) 'error in reading header of OUS file'
        stop 'error stop read_ous_header: reading header'
        end

c***************************************************************

	subroutine write_ous_info(iu)

	implicit none

	integer iu

	integer nvers,nkn,nel,nlv,nvar
	integer date,time
	real href,hzmin
	character*80 title,femver

        call getous(iu,nvers,nkn,nel,nlv)
        call ous_get_date(iu,date,time)
        call ous_get_title(iu,title)
        call ous_get_femver(iu,femver)
        call ous_get_hparams(iu,href,hzmin)

        write(6,*) 'nvers:  ',nvers
        write(6,*) 'nkn:    ',nkn
        write(6,*) 'nel:    ',nel
        write(6,*) 'nlv:    ',nlv
        !write(6,*) 'nvar:   ',nvar
        write(6,*) 'title:  ',trim(title)
        write(6,*) 'femver: ',trim(femver)
        write(6,*) 'date:   ',date
        write(6,*) 'time:   ',time
        write(6,*) 'href:   ',href
        write(6,*) 'hzmin:  ',hzmin

	end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine open_ous_type(type,status,nunit)

c open OUS file with default simulation name and given extension

        implicit none

        character*(*) type,status
        integer nunit

        integer nb
        character*80 file

        integer ifileo

        call def_make(type,file)
        nb = ifileo(0,file,'unform',status)

        if( nb .le. 0 ) then
          write(6,*) 'file: ',file
          stop 'error stop open_ous_type: opening file'
        end if

        nunit = nb

        end

c***************************************************************

        subroutine open_ous_file(name,status,nunit)

        implicit none

        character*(*) name,status
        integer nunit

        integer nb
        character*80 file

        integer ifileo

        call mkname(' ',name,'.ous',file)
        nb = ifileo(0,file,'unform',status)
        if( nb .le. 0 ) then
	  write(6,*) 'file: ',file
	  stop 'error stop open_ous_file: opening file'
	end if

        nunit = nb

        end

c***************************************************************

        subroutine qopen_ous_file(text,status,nunit)

c asks for name and opens ous file

        implicit none

        character*(*) text,status
        integer nunit

        character*80 name

        write(6,*) text
        read(5,'(a)') name
        write(6,*) name

        call open_ous_file(name,status,nunit)

        end

c***************************************************************
c***************************************************************
c***************************************************************

        subroutine ous_get_it_start(file,itstart)

c gets it of first record

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

c***************************************************************

        subroutine ous_get_it_end(file,itend)

c gets it of last record

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

c***************************************************************

        function check_ous_file(file)

        implicit none

        logical check_ous_file
        character*(*) file

        integer nb,nvers
        integer ifileo

        check_ous_file = .false.

        nb = ifileo(0,file,'unform','old')
        if( nb .le. 0 ) return
        call ous_is_ous_file(nb,nvers)
        close(nb)

        check_ous_file = nvers > 0
        
        end

c***************************************************************

        subroutine velzeta2scal(nel,nkn,nlv,nlvddi,nen3v,ilhkv
     +				,zenv,uprv,vprv
     +                          ,zv,sv,dv)

	implicit none

	integer nkn,nel,nlv,nlvddi
	integer nen3v(3,nel)
	integer ilhkv(nkn)
	real zenv(3,nel)
	real uprv(nlvddi,nkn)
	real vprv(nlvddi,nkn)
	real zv(nkn)
	real sv(nlvddi,nkn)
	real dv(nlvddi,nkn)

	integer ie,k,ii,l,lmax
	real u,v,s,d

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

c***************************************************************

