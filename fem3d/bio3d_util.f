c
c $Id: bio3d_util.f,v 1.3 2008-10-10 09:29:54 georg Exp $
c
c bio3d_util - utility routines for bio3d
c
c revision log :
c
c 18.04.2008    ggu     copied from weutro_sedim.f
c 09.10.2008    ggu     new call to confop
c 09.03.2012    ggu     bug fix: ilhkv was real
c
c********************************************************************

	subroutine loicz1(knode,vol,depth)

c EUTRO 0-D (LOICZ BUdgeting Procedure)
c
c new version -> does everything: initializes, accumulates, writes

	use basin

	implicit none

	integer nlzstate
	parameter(nlzstate=3)

        integer knode           !node for accumulation - 0 if init or write
        real vol                !volume [m**3]
        real depth              !depth of box [m]

        include 'param.h'
        include 'donata.h'


        real, save, allocatable :: elz(:,:)     !loicz budg proc ariables

        logical bloicz
        integer i,k
        real elzaux(nlzstate)    !diagnostic variable
        real tlztot(nlzstate)
        real getpar

        integer iublz,itmconlz,idtconlz
        save iublz,itmconlz,idtconlz

        integer icall
        save icall
        data icall / 0 /

c-----------------------------------------------------------
c see if routine has to be executed
c-----------------------------------------------------------

        bloicz = .true.

        if( icall .eq. -1 ) return

        if( .not. bloicz ) then
          icall = -1
          return
        end if

c-----------------------------------------------------------
c initialization
c-----------------------------------------------------------

        if( icall .eq. 0 ) then

	  allocate(elz(nkn,nlzstate))

          do i=1,nlzstate
            do k=1,nkn
              elz(k,i) = 0.
            end do
          end do

          iublz = 0
          itmconlz = nint(getpar('itmcon'))
          idtconlz = nint(getpar('idtcon'))

          call confop(iublz,itmconlz,idtconlz,1,nlzstate,'lcz')

          write(6,*) 'bio3d  loicz budget initialized...'

          icall = 1

          return
        end if 

c-----------------------------------------------------------
c accumulation of results
c-----------------------------------------------------------

        if( knode .gt. 0 ) then
          elzaux(1) = (prod - cons)*vol       !nem
          elzaux(2) = (ddin1+ddin2)*vol       !ddin
          elzaux(3) = denit*vol               !denitrificazione

          do i=1,nlzstate
            elz(knode,i) = elzaux(i)
          end do

          return
        end if

c-----------------------------------------------------------
c diagnostics and write
c-----------------------------------------------------------

        do i=1,nlzstate
          call scalmass(elz(1,i),depth,tlztot(i))   !mass ctrl loicz
        end do

        do i=1,nlzstate
          call confil(iublz,itmconlz,idtconlz,95+i,1,elz(1,i))
        end do

        call lcz_av_shell(elz)          !aver/min/max of nem and ddin

c-----------------------------------------------------------
c end of routine
c-----------------------------------------------------------

	end

c********************************************************************

	subroutine sed_av_shell(es)

c computes and writes average/min/max of sed variables
c
c id = 360
c
c es1) average	== 361
c es1) min	== 362
c es1) max	== 363
c es2) average	== 364
c ...

	use basin

	implicit none

c parameter

	include 'param.h'

	integer nsstate
	parameter( nsstate = 2 )

	real es(nkndi,nsstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: sedacu(:,:)
	real, save, allocatable :: sedmin(:,:)
	real, save, allocatable :: sedmax(:,:)

	integer ivect(8)
	save ivect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nsstate

	  allocate(sedacu(nkn,nsstate))
	  allocate(sedmin(nkn,nsstate))
	  allocate(sedmax(nkn,nsstate))

	  id = 360
	  call cmed_init('sdv',id,nvar,1,idtc,itmc
     +				,sedacu,sedmin,sedmax,ivect)

	  icall = 1
	end if

	call cmed_accum(1,es,sedacu,sedmin,sedmax,ivect)

	end

c********************************************************************

	subroutine lcz_av_shell(elz)

c computes and writes average/min/max of bio variables
c
c id = 460
c
c elz) average	== 461
c elz) min	== 462
c elz) max	== 463
c elz) average	== 464
c ...

	use basin

	implicit none

c parameter

	include 'param.h'

	integer nlzstate
	parameter( nlzstate = 3 )

	real elz(nkndi,nlzstate)	!state vector

c local
	integer idtc,itmc,itsmed
	integer id,nvar
c function
	real getpar
c save
	double precision, save, allocatable :: lczacu(:,:)
	real, save, allocatable :: lczmin(:,:)
	real, save, allocatable :: lczmax(:,:)

	integer ivect(8)
	save ivect

	integer icall
	save icall

	data icall / 0 /

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then

          itsmed=nint(getpar('itsmed'))
          if( itsmed .le. 0 ) then
            icall = -1
            return
          end if

	  idtc=nint(getpar('idtcon'))
	  itmc=nint(getpar('itmcon'))

	  nvar = nlzstate

	  allocate(lczacu(nkn,nlzstate))
	  allocate(lczmin(nkn,nlzstate))
	  allocate(lczmax(nkn,nlzstate))

	  write (6,*) 'cmed loicz inizializzato'
	  id = 460
	  call cmed_init('lzv',id,nvar,1,idtc,itmc
     +				,lczacu,lczmin,lczmax,ivect)

	  icall = 1
	end if

	!write(6,*) 'lzc: ',id,ivect(1)

	call cmed_accum(1,elz,lczacu,lczmin,lczmax,ivect)

	end

c********************************************************************

        subroutine setsedload(nlvddi,nknddi,nstate,eload,elini)

c sets up sediment loading

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        integer nlvddi,nknddi,nstate
        real eload(nlvddi,nknddi,nstate)
        real elini(nstate)

	include 'param.h'

        integer mode,i,k,l,lmax
        real d,vol,area

        mode = -1

        do i=1,nstate
         do k=1,nkn
          lmax = ilhkv(k)
          do l=1,lmax
            eload(l,k,i) = 0.
          end do
          call dvanode(lmax,k,mode,d,vol,area)   !gets depth, volume and area
          eload(lmax,k,i) = elini(i) * area
         end do
        end do

        end

c********************************************************************

        subroutine check_es(es)

	use basin

        include 'param.h'

        integer nsstate
        parameter( nsstate = 2 )

        real es(nkndi,nsstate)         !sediment state variables

        integer k,i

        i = 1

        k=801
        write(6,*) 'es: ',k,i,es(k,i)

        k=1201
        write(6,*) 'es: ',k,i,es(k,i)

        end

c********************************************************************

