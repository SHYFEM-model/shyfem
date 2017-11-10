c
c $Id: weutro_wdia.f,v 1.3 2015-11-25 donata e leslie $
c
c weutro_wdia - dianostic variables routines for weutro
c
c revision log :
c
c 10.10.2017    dmk&laa written from scratch
c 26.10.2017    ggu     integrated into main trunk
c
c********************************************************************

	subroutine wdia(k,t,dt,vol,depth,cdia)

c EUTRO 0-Dimensional (Sediments)

        implicit none

        integer nstate          !total number of state parameters
        parameter( nstate = 3 )

	integer k
	real t			!actual time [day]
        real dt                 !time step [day]
        real vol                !volume [m**3]
        real depth              !depth of box [m]
        real vel                !velocity [m/s]
c        real stp                !temperature [C]
        real cds(nstate)        !source term
        real cdia(nstate)          !output auxiliary variable  
        real cdiaold(nstate)       !old auxiliary variable
        real ca(nstate)            !auxiliary state vector     
        real caold (nstate)       !old auxiliary state vector

        real GPout,DPPout,GRZout

	logical debug
	logical bdia

	integer i
	real wpsink,wnsink,kpresusp,knresusp,kpcsed,kncsed,fvel
	real volsed,resusp
        real var1, var2, var3

        integer icall
        save icall
        data icall / 0 /

	bdia = .false.
	bdia = .true.

	if( .not. bdia ) return

        call leggiGPP(GPout)
        call leggiDPP(DPPout)
        call leggiGRZ(GRZout)

        var1=cdia(1)
        var2=cdia(2)
        var3=cdia(3)
        
        ca(1) = var1 !
        ca(2) = var2 !
        ca(3) = var3 !

         cds(1) = GPout*depth*(dt/86400)    !g C/m2/day gross primary production
c        cds(1) = GPout !output g C/m3/day

        cds(2) = (GPout-DPPout)*depth*(dt/86400)    !net primary production
        cds(3) = (dt/86400)   !DPPout*depth      !net primary production
        write(88,*) 'check:',ca,caold,cdia
cc
        ca(1) = cdia(1) +cds(1)
        ca(2) = cdia(2) +cds(2)
        ca(3) = cdia(3) +cds(3)  

        ! ca = (ca * depth + dt * cds)/depthnew

c       call euler(1,dt,vol,ca(1),caold,cds(1))
c       call euler(2,dt,vol,ca(2),caold(2),cds(2))
        write(88,*) 'check:',ca,caold,cdia

        cdia(1)=ca(1)
        cdia(2)=ca(2)
        cdia(3)=ca(3)

	end

c********************************************************************

        subroutine leggiGPP(gppread)

        implicit none
        include 'weutro.h'

        real gppread

        gppread=GPP

        return
        end

c********************************************************************

        subroutine leggiDPP(DPPread)

        implicit none
        include 'weutro.h'

        real DPPread

        DPPread=DPP

        return
        end

c********************************************************************

        subroutine leggiGRZ(GRZread)

        implicit none
        include 'weutro.h'

        real GRZread

        GRZread=GRZ

        return
        end

c********************************************************************

