c
c $Id: subwwm.f,v 1.4 2008-10-10 09:29:54 georg Exp $
c
c fifo pipe routines for WWM wave model
c
c routines for reading and writing data from the SHYFEM and the WWM model
c througth the FIFO PIPE mechanism
c file format: binary file with each value per node
c
c file unit to write:
c 	110 ---> velocity u				(velx)
c 	111 ---> velocity v				(vely)
c	112 ---> water level (z)			(lev)
c	113 ---> bathymetry (hkv) 			(bot)
c	114 ---> x wind component			(windx)
c	115 ---> y wind component			(windy)
c file unit to read:
c 	101 ---> x gradient of the radiation stresses	(stressx)
c 	102 ---> y gradient of the radiation stresses	(stressy)
c	103 ---> significant wave heigh			(waveh)
c	104 ---> significant wave period		(wavet)
c	105 ---> significant wave direction		(waved)
c	106 ---> orbital velocity			(orbit)
c	107 ---> stokes drift x				(stokesx)
c	108 ---> stokes drift y				(stokesy)
c
c revision log :
c
c 18.10.2006    ccf     integrated into main tree
c 10.04.2008    ccf     added waveov and stokesx/y
c 09.10.2008    ggu     new call to confop
c
c**************************************************************

        subroutine init_pipe(ipipe,idcoup)

c open pipes

	implicit none

        include 'param.h'

        integer ipipe		!call for coupling with pipe
	integer idcoup		!time step for sincronizing with wwm [s]

        real radx(neldim),rady(neldim)
        common /radx/radx,/rady/rady
        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave period [s]
        real waved(nkndim)      !wave direction
        real waveov(nkndim)	!orbital velocity
        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov

        real dtwave             !WWM wave model time step [s]
	integer iwave 		!call for wave model [2=wwm]
	real getpar		!get parameter function
	integer k,ie
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

c-------------------------------------------------------------
c initialize wave parameters
c-------------------------------------------------------------

        do ie = 1,neldim
          radx(ie) = 0.
          rady(ie) = 0.
        end do

	do k = 1,nkndim
          waveh(k) = 0.
          wavep(k) = 0.
          waved(k) = 0.
          waveov(k) = 0.
        end do

c-------------------------------------------------------------
c find out what to do
c-------------------------------------------------------------

	ipipe = 0
	idcoup = 0

        iwave = nint(getpar('iwave'))
        if( iwave .eq. 2 ) ipipe = 1
        if( ipipe .le. 0 ) return

c-------------------------------------------------------------
c open pipe files
c-------------------------------------------------------------
        
        open(110,file='p_velx.dat',form='unformatted')
        open(111,file='p_vely.dat',form='unformatted')
        open(112,file='p_lev.dat',form='unformatted')
        open(113,file='p_bot.dat',form='unformatted')
        open(114,file='p_windx.dat',form='unformatted')
        open(115,file='p_windy.dat',form='unformatted')

        open(101,file='p_stressx.dat',form='unformatted')
        open(102,file='p_stressy.dat',form='unformatted')
        open(103,file='p_waveh.dat',form='unformatted')
        open(104,file='p_wavet.dat',form='unformatted')
        open(105,file='p_waved.dat',form='unformatted')
        open(106,file='p_orbit.dat',form='unformatted')
        open(107,file='p_stokesx.dat',form='unformatted')
        open(108,file='p_stokesy.dat',form='unformatted')

c-------------------------------------------------------------
c set coupling time step 
c-------------------------------------------------------------
       
        dtwave = getpar('dtwave')
	idcoup = dtwave

        write(6,*) 'spectral wave module has been initialized'

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------
        
	end

c**************************************************************

	subroutine read_pipe(ipipe,it,idcoup)

c read from PIPE

        implicit none

        include 'param.h'	

        integer ipipe
        integer it		!time
	integer idcoup		!time step for sincronizing with wwm [s]

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer k

        real waveh(nkndim)      !wave height [m]
        real wavep(nkndim)      !wave period [s]
        real waved(nkndim)      !wave direction
        real waveov(nkndim)	!orbital velocity
        real stokesx(nkndim)	!stokes velocity x
        real stokesy(nkndim)	!stokes velocity y
        common /waveh/waveh, /wavep/wavep, /waved/waved, /waveov/waveov
	common /stokesx/stokesx, /stokesy/stokesy

	real forcex(nkndim)	!gradient of the radiation stress
        real forcey(nkndim)

        real radx(neldim),rady(neldim)
        common /radx/radx,/rady/rady

	real getpar		!get parameter function
        integer ius,itmcon,idtcon
        save ius,itmcon,idtcon
 
        integer icall           !initialization parameter                        
        save icall
        data icall /0/

        if (ipipe .le. 0 ) return

c       -----------------------------------------------
c       open output file for waves
c       -----------------------------------------------

        if( icall .eq. 0 ) then
            itmcon = nint(getpar('itmcon'))
            idtcon = nint(getpar('idtcon'))
            call confop(0,itmcon,idtcon,1,3,'wav')
            icall = 1
        end if
                    
c       -----------------------------------------------
c       same time step, fix me
c       -----------------------------------------------

        if (mod(it,idcoup) .eq. 0 ) then
 
c           -----------------------------------------------
c           read stress and wave characteristics
c           -----------------------------------------------

            do k = 1,nkn 
              read(101) forcex(k)
              read(102) forcey(k)
              read(103) waveh(k)
              read(104) wavep(k)
              read(105) waved(k)
              read(106) waveov(k)
              read(107) stokesx(k)
              read(108) stokesx(k)
            end do

c           -----------------------------------------------
c           writes output to the file.wav 
c           -----------------------------------------------
            
            call confil(ius,itmcon,idtcon,31,1,waveh)
            call confil(ius,itmcon,idtcon,32,1,wavep)
            call confil(ius,itmcon,idtcon,33,1,waved)
            
c           -----------------------------------------------
c           convert node stress value to element value
c           -----------------------------------------------

            call n2e2d(forcex,radx)
            call n2e2d(forcey,rady)

            write(*,*) 'SHYFEM read stress and wave ',it

        end if

        end

c**************************************************************

        subroutine write_pipe(ipipe,it,idcoup)

c write to PIPE

        implicit none

        include 'param.h'

        integer ipipe
        integer it              !time [s]
	integer idcoup		!time step for sincronizing with wwm [s]

        integer k
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        real wxnv(nkndim),wynv(nkndim)    !x and y wind component [m/s]
        common /wxnv/wxnv,/wynv/wynv
        real uprv(nlvdim,1), vprv(nlvdim,1)
        common /uprv/uprv, /vprv/vprv
        real znv(1)
        common /znv/znv
        real hkv(1)
        common /hkv/hkv

        if (ipipe .le. 0 ) return

        if (mod(it,idcoup) .eq. 0 ) then

c           -----------------------------------------------
c           write velocities and water level
c           -----------------------------------------------

            do k = 1,nkn 
              write(110) uprv(1,k)
              write(111) vprv(1,k)
              write(112) znv(k)
              write(113) hkv(k)
              write(114) wxnv(k)
              write(115) wynv(k)
            end do

            write(*,*) 'SHYFEM write vel and water level', it
 
        end if

        end

c**************************************************************************

c first is HS, Period, Wave Direction, Directional Spreading = cos^MS you give MS 
c and the spectral type is ...
c 1 PM spectrum
c 2 JONSWAP spectrum
c 3 Single bin
c 4 Gaussian Distribution

