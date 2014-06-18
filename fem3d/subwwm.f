!
! $Id: subwwm.f,v 1.1 2006/10/18 15:35:13 georg Exp $
!
! fifo pipe routines for WWM wave model
!
! routines for reading and writing data from the SHYFEM and the WWM model
! througth the FIFO PIPE mechanism
! file format: binary file with each value per node
!
! IWAVE:
!   - 1: the parametric wave model is called (see subwave.f)
!   - >2: the spectral wave model is called (WWM)
!        - 2: wind is elaborated by SHYFEM and passed to WWM 
!	 - 3: wind is elaborated by WWM and passed to SHYFEM 
!
! file unit to write:
! 	120 ---> velocity u
! 	121 ---> velocity v
!	122 ---> water level (z)
!	123 ---> bathymetry (hkv) & number of layers
!	124 ---> x wind component
!	125 ---> y wind component
!       126 ---> 3D layer depth (in the middle of layer)
! file unit to read:
! 	101 ---> x gradient of the radiation stresses (radx)
! 	102 ---> y gradient of the radiation stresses (rady)
!	103 ---> significant wave heigh
!	104 ---> mean period
!	105 ---> significant wave direction
!	106 ---> KME
!	107 ---> peak period
!	108 ---> KP
!	109 ---> orbital velocity
!	110 ---> stokes drift x
!	111 ---> stokes drift y
!
! revision log :
!
! 18.10.2006    ccf     integrated into main tree
! 19.06.2008    aac&ccf udate to 3D version
! 16.04.2009	ccf	update to new WWMII-2 version, both 2D and 3D
! 18.02.2011	ggu	compiler warnings/errors adjusted
!
!**************************************************************

        subroutine init_pipe(ipipe,idcoup)

! open pipes

	implicit none

        include 'param.h'

        integer ipipe		!call for coupling with pipe
	integer idcoup		!time step for sincronizing with wwm [s]

        real radx(nlvdim,neldim),rady(nlvdim,neldim)
        common /radx/radx,/rady/rady
        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave period [s]
        real waved(nkndim)      !wave direction
        real waveov(nkndim)	!orbital velocity
        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov

	real ddl(nlvdim,nkndim)
	common /ddl/ddl

        integer ius,itmcon,idtcon
	integer iwave 		!call for wave model [2=wwm]
	real getpar		!get parameter function
	integer k,ie,l

!-------------------------------------------------------------
! initialize wave parameters
!-------------------------------------------------------------

        do ie = 1,neldim
  	  do l = 1,nlvdim
            radx(l,ie) = 0.
            rady(l,ie) = 0.
	  end do
        end do

	do k = 1,nkndim
          waveh(k) = 0.
          wavep(k) = 0.
          waved(k) = 0.
          waveov(k) = 0.
	  do l = 1,nlvdim
	    ddl(l,k) = 0.
	  end do
        end do

!-------------------------------------------------------------
! find out what to do
!-------------------------------------------------------------

	idcoup = 0

        iwave = nint(getpar('iwave'))
        if( iwave .eq. 2 ) then
	  ipipe = 1	!wind from SHYFEM
	elseif ( iwave .eq. 3) then
	  ipipe = 2	!wind from WWM
	else
	  ipipe = 0	!no SHYFEM-WWM coupling
	end if

        if( ipipe .le. 0 ) return

!-------------------------------------------------------------
! open pipe files
!-------------------------------------------------------------
        
        open(120,file='p_velx.dat',form='unformatted')
        open(121,file='p_vely.dat',form='unformatted')
        open(122,file='p_lev.dat',form='unformatted')
        open(123,file='p_bot.dat',form='unformatted')
        open(124,file='p_windx.dat',form='unformatted')
        open(125,file='p_windy.dat',form='unformatted')
        open(126,file='p_zeta3d.dat',form='unformatted')

        open(101,file='p_stressx.dat',form='unformatted')
        open(102,file='p_stressy.dat',form='unformatted')
        open(142,file='p_stresxy.dat',form='unformatted')
        open(103,file='p_waveh.dat',form='unformatted')
        open(104,file='p_wavet.dat',form='unformatted')
        open(105,file='p_waved.dat',form='unformatted')
        open(106,file='p_wavekm.dat',form='unformatted')
        open(107,file='p_wavetp.dat',form='unformatted')
        open(108,file='p_wavekp.dat',form='unformatted')
        open(109,file='p_orbit.dat',form='unformatted')
        open(110,file='p_stokesx.dat',form='unformatted')
        open(111,file='p_stokesy.dat',form='unformatted')

!-------------------------------------------------------------
! open output file 
!-------------------------------------------------------------

        ius = 31
        itmcon = nint(getpar('itmcon'))
        idtcon = nint(getpar('idtcon'))
        call confop(ius,itmcon,idtcon,1,3,'wav')

!-------------------------------------------------------------
! set coupling time step 
!-------------------------------------------------------------

        idcoup = nint(getpar('dtwave'))

        write(6,*) 'spectral wave module has been initialized'
 	!call getwwmbound

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------
        
	end

!**************************************************************

	subroutine read_pipe(ipipe,it,idcoup)

! read from PIPE

        implicit none

        include 'param.h'	

        integer ipipe
        integer it		!time
	integer idcoup		!time step for sincronizing with wwm [s]

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer k,l

        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave mean period [s]
        real waved(nkndim)      !wave direction
        real wavekm(nkndim)     !wave KME
        real wavepp(nkndim)     !wave peak period [s]
        real wavekp(nkndim)     !wave KP
        real waveov(nkndim)	!wave orbital velocity
        real stokesx(nkndim)	!stokes velocity x
        real stokesy(nkndim)	!stokes velocity y
        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov
	common /stokesx/stokesx, /stokesy/stokesy

	real forcex(nlvdim,nkndim)	!gradient of the radiation stress
        real forcey(nlvdim,nkndim)
	real SXX3D(nlvdim,nkndim)
	real SYY3D(nlvdim,nkndim)
	real SXY3D(nlvdim,nkndim)

        real radx(nlvdim,neldim),rady(nlvdim,neldim)
        common /radx/radx,/rady/rady

        real wxv(nkndim),wyv(nkndim)    !x and y wind component [m/s]
        common /wxv/wxv,/wyv/wyv
        real tauxnv(nkndim),tauynv(nkndim)
        common /tauxnv/tauxnv,/tauynv/tauynv
        real ppv(nkndim)
        common /ppv/ppv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

	real getpar		!get parameter function
        integer ius,itmcon,idtcon
        save ius,itmcon,idtcon

        integer icall           !initialization parameter                        
        save icall
        data icall /0/

        if (ipipe .le. 0 ) return

!       -----------------------------------------------
!       open output file for waves
!       -----------------------------------------------

        if( icall .eq. 0 ) then
            ius = 31
            itmcon = nint(getpar('itmcon'))
            idtcon = nint(getpar('idtcon'))
            do k = 1,nkndim
              do l = 1,nlvdim
                 forcex(l,k) = 0.
                 forcey(l,k) = 0.
		 SXX3D(l,k) = 0.
		 SXY3D(l,k) = 0.
		 SYY3D(l,k) = 0.
              end do
            end do
            icall = 1
        end if

!       -----------------------------------------------
!       same time step, fix me
!       -----------------------------------------------

        if (mod(it,idcoup) .eq. 0 ) then
 
!         -----------------------------------------------
!         read stress and wave characteristics
!         -----------------------------------------------

  	  if (ipipe .eq. 1) then 	!do not read wind from WWM

            do k = 1,nkn

	      do l = 1,nlv 
                 !read(101) forcex(l,k)
                 !read(102) forcey(l,k)
                 read(101) SXX3D(l,k)
                 read(102) SYY3D(l,k)
                 read(142) SXY3D(l,k)
	      end do

              read(103) waveh(k)
              read(104) wavep(k)
              read(105) waved(k)
              read(106) wavekm(k)
              read(107) wavepp(k)
              read(108) wavekp(k)
              read(109) waveov(k)
              read(110) stokesx(k)
              read(111) stokesy(k)

            end do

	  elseif (ipipe .eq. 2) then	!read wind from WWM

            do k = 1,nkn

	      do l = 1,nlv 
                 !read(101) forcex(l,k)
                 !read(102) forcey(l,k)
                 read(101) SXX3D(l,k)
                 read(102) SYY3D(l,k)
                 read(142) SXY3D(l,k)
	      end do

              read(103) waveh(k)
              read(104) wavep(k)
              read(105) waved(k)
              read(106) wavekm(k)
              read(107) wavepp(k)
              read(108) wavekp(k)
              read(109) waveov(k)
              read(110) stokesx(k)
              read(111) stokesy(k)
              read(124) wxv(k)
              read(125) wyv(k)
              ppv(k) = 0.

            end do

            call wstress(nkn,wxv,wyv,tauxnv,tauynv)

	  end if

!         -----------------------------------------------
!         check for NaN in radiation stresses
!         -----------------------------------------------

	  !call nantest(nkn*nlvdim,forcex,'forcex')
	  !call nantest(nkn*nlvdim,forcey,'forcey')
	  call nantest(nkn*nlvdim,SXX3D,'SXX3D')
	  call nantest(nkn*nlvdim,SXY3D,'SXY3D')
	  call nantest(nkn*nlvdim,SYY3D,'SYY3D')

!         -----------------------------------------------
!         convert node stress value to element value
!         -----------------------------------------------
 
          !call n2e3d(nlvdim,forcex,radx)
          !call n2e3d(nlvdim,forcey,rady)

	  call diffxy(SXX3D,SYY3D,SXY3D,radx,rady)

          write(*,*) 'SHYFEM read stress and wave ',it

        end if

!       -----------------------------------------------
!       writes output to the file.wav 
!       -----------------------------------------------
            
        call confil(ius,itmcon,idtcon,31,1,waveh)
        call confil(ius,itmcon,idtcon,32,1,wavep)
        call confil(ius,itmcon,idtcon,33,1,waved)
 
        end

!**************************************************************

        subroutine write_pipe(ipipe,it,idcoup)

! write to PIPE

        implicit none

        include 'param.h'

        integer ipipe
        integer it              !time [s]
	integer idcoup		!time step for sincronizing with wwm [s]

        integer k,l,nlev,lmax
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real wxv(nkndim),wyv(nkndim)    !x and y wind component [m/s]
        common /wxv/wxv,/wyv/wyv
        real up0v(nkndim), vp0v(nkndim)
        common /up0v/up0v, /vp0v/vp0v
        real znv(nkndim)
        common /znv/znv
        real hkv(nkndim)
        common /hkv/hkv
        integer ilhkv(nkndim)           !number of node level
        common /ilhkv/ilhkv
        integer nlvdi,nlv
        common /level/ nlvdi,nlv
	real ddl(nlvdim,nkndim)		!3D layer depth (in the middle of layer)
	common /ddl/ddl
	real h(nlvdim)

        if (ipipe .le. 0 ) return

        if (mod(it,idcoup) .eq. 0 ) then

          do k = 1,nkn
	    call dep3dnod(k,+1,nlev,h)
            ddl(1,k) = - 0.5 * h(1)
            do l = 2,nlev
              ddl(l,k) = ddl(l-1,k) - 0.5 * (h(l) + h(l-1))
            end do
          end do

!         -----------------------------------------------
!         write velocities and water level
!         -----------------------------------------------

	  if (ipipe .eq. 1) then	! write wind to wwm

            do k = 1,nkn 
              write(120) up0v(k)
              write(121) vp0v(k)
              write(122) znv(k)
              write(123) hkv(k)
              write(123) ilhkv(k)
              write(124) wxv(k)
              write(125) wyv(k)

              do l = 1,nlv
                write(126) ddl(l,k)
	      end do 
	    end do 

	  elseif (ipipe .eq. 2) then	! do not write wind to wwm

            do k = 1,nkn 
              write(120) up0v(k)
              write(121) vp0v(k)
              write(122) znv(k)
              write(123) hkv(k)
              write(123) ilhkv(k)

              do l = 1,nlv
                write(126) ddl(l,k)
	      end do 
            end do

	  end if

          write(*,*) 'SHYFEM write vel and water level', it
 
        end if

        end

!**************************************************************************

        subroutine getwwmbound

!routine to write boundary node file for wwm (fort.22)

        implicit none

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv
        real hkv(1)
        common /hkv/hkv
	integer ibc,k,knode,kranf,krend,nn

	integer kbnd

        do ibc=1,nbc
	  call kanfend(ibc,kranf,krend)
	  nn = krend-kranf+1

	  write(26,*)nn

          do k=kranf,krend
            !knode=irv(k)
            knode=kbnd(k)
            write(26,25)knode,xgv(knode),ygv(knode),hkv(knode)
	  end do
	enddo

25	format(i10,3e14.4)

	end

!**************************************************************************
!differenzation of radiation stresses

        subroutine diffxy(SXX3D,SYY3D,SXY3D,radx,rady)

        implicit none

        include 'param.h'
        include 'evmain.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhv(neldim)
        common /ilhv/ilhv
        integer nen3v(3,neldim)        !node number
        common /nen3v/nen3v
        real radx(nlvdim,neldim),rady(nlvdim,neldim)
        double precision b,c           !x and y derivated form function [1/m]
	integer k,ie,ii,l,ilevel
	real radsx,radsy

	real SXX3D(nlvdim,nkndim)
	real SYY3D(nlvdim,nkndim)
	real SXY3D(nlvdim,nkndim)

	do ie = 1,nel
	  ilevel = ilhv(ie)
	  do l=1,ilevel
	    radsx = 0.
	    radsy = 0.
	    do ii = 1,3
	      k = nen3v(ii,ie)
              b = ev(3+ii,ie)
              c = ev(6+ii,ie)
	      radsx = radsx -(SXX3D(l,k)*b + SXY3D(l,k)*c)
	      radsy = radsy -(SXY3D(l,k)*b + SYY3D(l,k)*c)
	    end do
	  radx(l,ie) = radsx
	  rady(l,ie) = radsy
          end do
	enddo

	end

!**************************************************************************
