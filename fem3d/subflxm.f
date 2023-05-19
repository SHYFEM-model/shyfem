
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2002-2003,2006-2007,2006-2007  Georg Umgiesser
!    Copyright (C) 2009-2020  Georg Umgiesser
!    Copyright (C) 2017  Christian Ferrarin
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

c subroutines for flux computation - mpi versions
c
c revision log :
c
c 28.05.2022	ggu	copied from subflxa.f
c 19.01.2023	ggu	new routine gather_sum_d3()
c 21.04.2023	ggu	problems in gather_sum_d3() ... not yet solved
c 18.05.2023	ggu	flux_write() changed, new routine flx_collect_3d()
c
c******************************************************************
c******************************************************************
c******************************************************************

	subroutine correct_nlayers(nsect,nlayers,nlmax)

! corrects nlayers between domains

	use shympi

	implicit none

	integer nsect			!total number of sections
	integer nlayers(nsect)		!layers of sections
	integer nlmax			!maximum layers in all sections

	integer is
	!integer, allocatable :: nlayers_domain(:,:)
	integer nlayers_domain(nsect,n_threads)

	!allocate(nlayers_domain(nsect,n_threads))

	call shympi_gather(nlayers,nlayers_domain)

	do is=1,nsect
	  nlayers(is) = maxval(nlayers_domain(is,:))
	end do

	nlmax = maxval(nlayers)

	end

c******************************************************************

	subroutine correct_iflux(kfluxm,kflux,iflux)

	use shympi

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)

	integer i,k,kext
	integer ipext

	do i=1,kfluxm
	  k = kflux(i)
	  if( k <= 0 ) cycle
	  if( id_node(k) /= my_id ) then
	    kext = ipext(k)
	    write(6,*) 'deleting section node: ',my_id,id_node(k),kext
	    iflux(1,i) = 0
	  end if
	end do

	end

c******************************************************************

        subroutine convert_nodes(n,nodes)

c converts external to internal nodes

        implicit none

	integer n
	integer nodes(n)

        integer i
        integer kext,kint
        integer ifirst,ilast,nnode

        integer ipint
        logical nextline

        if( n <= 0 ) return

        nnode = 0

        do while( nextline(nodes,n,nnode,ifirst,ilast) )
          do i=ifirst,ilast
            kext = nodes(i)
            kint = ipint(kext)
            if( kint > 0 ) then
              nodes(i) = kint
            else
              nodes(i) = -1
            end if
          end do
        end do

        end

c******************************************************************
c******************************************************************
c******************************************************************

	subroutine flux_initialize(kfluxm,kflux,iflux
     +				,nsect,nlayers,nlmax)

	use shympi

	implicit none

	integer kfluxm
	integer kflux(kfluxm)
	integer iflux(3,kfluxm)
	integer nsect
	integer nlayers(nsect)
	integer nlmax

	integer nsaux		! not needed

	call flx_init(kfluxm,kflux,nsaux,iflux)	!basic checks and iflux
	call shympi_barrier
	call get_nlayers(kfluxm,kflux,nsect,nlayers,nlmax)

	! now correct for mpi (exchange betrween domains)

	call shympi_barrier
	call correct_nlayers(nsect,nlayers,nlmax)
	call shympi_barrier
	call correct_iflux(kfluxm,kflux,iflux)

	end

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flux_file_open(type,da_out,nvar,nsect,kfluxm
     +				,kflux_ext,nlayers,chflx)

	use shympi

	implicit none

	character*(*) type
	double precision da_out(4)
	integer nvar
	integer nsect
	integer kfluxm
	integer kflux_ext(kfluxm)
	integer nlayers(nsect)
	character*80 chflx(nsect)

	include 'simul.h'

	logical debug
	integer nbflx,nvers,idtflx,nlmax
	integer ierr
	double precision atime0
	character*80 title,femver

	integer ifemop

	debug = .false.

	if( shympi_is_master() ) then

          nbflx=ifemop(type,'unform','new')
          if(nbflx.le.0) goto 99
	  da_out(4) = nbflx
	  nlmax = maxval(nlayers)

	  nvers = 0
	  idtflx = nint(da_out(1))
	  call flx_write_header(nbflx,nvers,nsect,kfluxm,idtflx
     +                                  ,nlmax,nvar,ierr)
	  if( ierr /= 0 ) goto 98

	  title = descrp
          call get_shyfem_version_and_commit(femver)
          call get_absolute_ref_time(atime0)

	if( debug ) then
	write(6,*) 'ggguuu3:'	!GGUFLUX
	write(6,*) nbflx,nvers,nsect,kfluxm
	write(6,*) nlayers
	write(6,*) atime0
	write(6,*) title
	write(6,*) femver
	write(6,*) chflx
	write(6,*) 'calling header2:'
	end if

          call flx_write_header2(nbflx,nvers,nsect,kfluxm
     +                          ,kflux_ext,nlayers
     +                          ,atime0,title,femver,chflx,ierr)
	  if( ierr /= 0 ) goto 98

	end if

	call shympi_barrier

	return
   99   continue
        write(6,*) 'Error opening FLX file :'
        stop 'error stop flux_file_open: opening flx file'
   98   continue
        write(6,*) 'Error writing header of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop flux_file_open: writing flx header'
	end

!******************************************************************

	subroutine flux_write(nbflx,atime,ivar,nlvddi,nsect
     +				,nlayers,flux_global)

	use shympi

	implicit none

	integer nbflx			!unit number
	double precision atime
	integer ivar
	integer nlvddi
	integer nsect
	integer nlayers(nsect)
	double precision flux_global(0:nlv_global,3,nsect)

	integer nvers
	integer ierr
	real flux_real(0:nlv_global,3,nsect)

	if( nlvddi /= nlv_global ) stop 'error stop flx_write: (1)'

	flux_real = flux_global

	nvers = 0

	if( shympi_is_master() ) then
	  call flx_write_record(nbflx,nvers,atime,nlv_global,nsect,ivar
     +				,nlayers,flux_real,ierr)
	  if( ierr /= 0 ) goto 97
	end if

	call shympi_barrier

	return
   97   continue
        write(6,*) 'Error writing file FLX'
        write(6,*) 'unit,ierr,ivar :',nbflx,ierr,ivar
        stop 'error stop flux_write: writing flx record'
	end

!******************************************************************

	subroutine flx_collect_2d(n,flux2d)

	use shympi

	implicit none

	integer n
	double precision flux2d(n)

	double precision flux_domain(n,n_threads)

	call shympi_gather(flux2d,flux_domain)
	flux2d(:) = SUM(flux_domain,dim=2)
	
	end

!******************************************************************

	subroutine flx_collect_3d(nlvddi,n,flux3d,flux3dg)

	use shympi

	implicit none

	integer nlvddi
	integer n
	double precision, intent(in) :: flux3d(0:nlvddi,n)
	double precision, intent(out) :: flux3dg(0:nlv_global,n)

	double precision flux_domain(0:nlv_global,n,n_threads)
	double precision flux_aux(0:nlv_global,n)

	flux_aux = 0.
	flux_aux(0:nlvddi,:) = flux3d
	call shympi_gather(flux_aux,flux_domain)
	flux3dg(:,:) = SUM(flux_domain,dim=3)
	
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine reshape_local(ns,nt,n,source,target)

	implicit none

	integer ns
	integer nt
	integer n
	double precision source(ns,n)
	double precision target(nt,n)

	integer i

	if( ns > nt ) stop 'error stop reshape_local: ns>nt'

	target = 0.

	do i=1,n
	  target(1:ns,i) = source(1:ns,i)
	end do

	end

!******************************************************************

	subroutine gather_sum_dr(ns,nt,n,source,target)

	use shympi

	implicit none

	integer ns,nt,n
	double precision source(ns,n)
	real target(nt*n)

	double precision flux_lin(nt*n)
	double precision flux_domain(nt*n,n_threads)

	call reshape_local(ns,nt,n,source,flux_lin)
	call shympi_gather(flux_lin,flux_domain)
	target(:) = SUM(flux_domain,dim=2)

	end

!******************************************************************

	subroutine gather_sum_dd(ns,nt,n,source,target)

	use shympi

	implicit none

	integer ns,nt,n
	double precision source(ns,n)
	double precision target(nt*n)

	double precision flux_lin(nt*n)
	double precision flux_domain(nt*n,n_threads)

	call reshape_local(ns,nt,n,source,flux_lin)
	call shympi_gather(flux_lin,flux_domain)
	target(:) = SUM(flux_domain,dim=2)

	end

!******************************************************************

	subroutine gather_max_i(n,nodes)

	use shympi

	implicit none

	integer n
	integer nodes(n)

	integer nodes_domain(n,n_threads)

	if( .not. bmpi ) return

	call shympi_gather(nodes,nodes_domain)
	nodes(:) = maxval(nodes_domain,dim=2)

	end

!******************************************************************

	subroutine gather_max_r(n,nodes)

	use shympi

	implicit none

	integer n
	real nodes(n)

	real nodes_domain(n,n_threads)

	if( .not. bmpi ) return

	call shympi_gather(nodes,nodes_domain)
	nodes(:) = maxval(nodes_domain,dim=2)

	end

!******************************************************************

	subroutine gather_sum_d2(n,vals)

	use shympi

	implicit none

	integer n
	double precision vals(n)

	double precision vals_domain(n,n_threads)

	if( .not. bmpi ) return

	call shympi_gather(vals,vals_domain)
	vals(:) = SUM(vals_domain,dim=2)

	end

!******************************************************************

	subroutine gather_sum_d3(nlvddi,n,vals)

	use shympi

	implicit none

	integer nlvddi
	integer n
	double precision vals(0:nlvddi,n)

	integer iu,l,i
	!double precision vals_domain(nlvddi,n,n_threads)
	double precision val_aux(0:nlv_global,n)
	double precision vals_domain(0:nlv_global,n,n_threads)

	if( .not. bmpi ) return

	vals_domain = 0.
	call shympi_gather(vals,vals_domain)
	
	val_aux(:,:) = SUM(vals_domain,dim=3)
	vals(0:nlvddi,:) = val_aux(0:nlvddi,:)

	end

!******************************************************************

