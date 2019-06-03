
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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

c weutro_wdia - dianostic variables routines for weutro
c
c revision log :
c
c 10.10.2017	dmk&laa	written from scratch
c 26.10.2017	ggu	integrated into main trunk
c 04.11.2017	ggu	changed VERS_7_5_34
c 03.04.2018	ggu	changed VERS_7_5_43
c 16.02.2019	ggu	changed VERS_7_5_60
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

        write(*,*) 'wdia entrato'
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
        write(*,*) 'wdia check1:',ca,caold,cdia
        ca(1) = cdia(1) +cds(1)
        ca(2) = cdia(2) +cds(2)
        ca(3) = cdia(3) +cds(3)  

        ! ca = (ca * depth + dt * cds)/depthnew

c       call euler(1,dt,vol,ca(1),caold,cds(1))
c       call euler(2,dt,vol,ca(2),caold(2),cds(2))
        write(88,*) 'check:',ca,caold,cdia

        write(*,*) 'wdia check2:',ca,caold,cdia
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

	subroutine weutro_check(text)

	use eutro

	implicit none

	character*(*) text

	write(6,*) trim(text)
	write(6,*) '  pelagic 1 2 3 4 7 8'
        write(*,*) SUM(e(:,:,1)),SUM(e(:,:,2)),
     +                   SUM(e(:,:,3)),SUM(e(:,:,4)),
     +                   SUM(e(:,:,7)),SUM(e(:,:,8)),
	write(6,*) '  ulva 1 2'
	write(6,*) SUM(eul(:,1)),SUM(eul(:,2))

	end

c********************************************************************

