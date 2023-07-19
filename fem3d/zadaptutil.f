
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2017,2019  Georg Umgiesser
!    Copyright (C) 2022  Luca Arpaia
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

c routines for handling z-adaptive layers
c
c revision log :
c
c 05.06.2023    lrp     introduce z-star
c 18.07.2023	lrp	rzmov read from shy
c
c notes:
c this file is used also in:      
c      
c	compute_zadaptive_info (subele.f)
c	get_zadapt_info (newexpl.f)
c	set_zadapt_info (subele.f)
c	get_zadaptivelayer_thickness (subele.f)
c	init_zadaptation (shyfem.f)
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module zadapt
!==================================================================

	implicit none

        integer, save, allocatable :: nadapt_com(:,:) !number of adaptive layers
        real   , save, allocatable :: hadapt_com(:,:) !closing depth of adaptive layers
	real   , save :: rzmov_com		      !parameter for moving surface layers

!==================================================================
        end module zadapt
!==================================================================

c******************************************************************

        subroutine get_rzmov_info(rzmov)

        use zadapt

        implicit none

        real rzmov

        rzmov = rzmov_com

        end

c******************************************************************

        subroutine set_rzmov_info(rzmov)

        use zadapt

        implicit none

        real rzmov

        rzmov_com = rzmov

        end

c******************************************************************

        subroutine get_zadapt_info(ie,nadapt,hadapt)

        use zadapt

        implicit none

	integer ie
        integer nadapt(4)
        real hadapt(4)

	!call check_sigma_initialized  !lrp do some check here

        nadapt = nadapt_com(:,ie)
        hadapt = hadapt_com(:,ie)

        end

c******************************************************************

        subroutine set_zadapt_info(ie,nadapt,hadapt)

        use zadapt

        implicit none

        integer ie
        integer nadapt(4)
        real hadapt(4)

        nadapt_com(:,ie) = nadapt
        hadapt_com(:,ie) = hadapt

        end

c******************************************************************

	subroutine compute_zadapt_info(z,hlv,nsig,lmax,lmin,nzad,hzad)

c returns relevant info for z-surface-adaptive layers

	implicit none

	real z			!water level
	real hlv(lmax)		!layer structure
	integer nsig		!total number of sigma layers
	integer lmax            !bottom layer  index
	integer lmin            !surface layer  index (return)
	integer nzad		!number of z-adaptive layers (return)
	real hzad		!closing depth of z-adaptive layers (return)

	integer l,levmax
	real getpar,rgridtop,rgridmov

        lmin = 1
        levmax = 0       !no adapation -> all to zero	

	if( nsig .eq. 0 ) then   

c---------------------------------------------------------
c surface layer index
c---------------------------------------------------------  

c         !for now commented: lmin = 1

c         rgridtop = getpar('rztop')
c      	  do l=1,lmax              !a threshold is used
c            if((-hlv(l)+rgridtop*(hlv(l)-hlv(l-1))).le.z) exit
c          end do
c          lmin=min(l,lmax)         !safety min: jlhv>=ilhv    

c---------------------------------------------------------
c lowest index of adaptive deforming layers
c--------------------------------------------------------- 

          call get_rzmov_info(rgridmov) !getpar('rzmov')
          do l=lmin,lmax-1 !-1 to skip bottom layer
            if(z.le.(-hlv(l)+rgridmov*(hlv(l)-hlv(l-1)))) then
              levmax = l+1
            else
              exit
            end if
          end do

	end if 

c---------------------------------------------------------
c compute nzad, hzad
c--------------------------------------------------------- 

        hzad = hlv(levmax)
        nzad = max(0,levmax-lmin+1) !+1 (min 2 adaptive layer)

	end

c******************************************************************

        subroutine compute_zadaptive_info(ie,nlv,lmin,lmax,hlv,z,htot,
     +					  nadapt,ladapt,hadapt,hdl)

c returns z-surface-adaptive layers info. Info is computed by node of element:
c number of adaptive layers       by node (3) + by ele (1)
c lowest index of adaptive layer  by node (3) + by ele (1)
c closing depth of adaptive layer by node (3) + by ele (1)
c coefficients of adaptive layers 

        implicit none

        integer ie              !element index
        integer nlv             !total number of layers
        integer lmin(3)         !top layer index	
	integer lmax		!bottom layer index
	real hlv(nlv)           !layer structure
        real z(3)               !water level
	real htot(3)		!water depth	
        integer nadapt(4)       !total number of adaptive layers (return)
	integer ladapt(4)       !lowest index of adaptive layers (return)
        real hadapt(4)          !closing depth of adaptive layers(return)
	real hdl(nlv,3)         !coefficient (return)

	integer l,ii,levmax,lmine,nsigma,nlev
	real getpar,r,hsigma,htop,hbot
	real den,check
        logical bsigma

	nadapt = 0         !no adaptation -> all to zero
	ladapt = 0	   !no adaptation -> all to zero
	hadapt = 0	   !no adaptation -> all to zero	
	hdl = 0.           !no adaptation -> all to zero

	r = getpar('rzmov')

        call get_sigma_info(nlev,nsigma,hsigma)
        bsigma = nsigma .gt. 0

	lmine = maxval(lmin)

        if( .not.bsigma ) then

c---------------------------------------------------------
c loop over nodes: adaptation is node-driven
c---------------------------------------------------------

	do ii=1,3
	  levmax = 0       !no adapation -> all to zero
	  do l=lmin(ii),lmax-1 !-1 to skip bottom layer

c---------------------------------------------------------
c lowest index of adaptive deforming layers
c---------------------------------------------------------	  
	  
            if(z(ii).le.(-hlv(l)+r*(hlv(l)-hlv(l-1)))) then 
              levmax = l+1

c---------------------------------------------------------
c no deformation
c--------------------------------------------------------- 

	    else  
	      exit 
	    end if  
          end do

c---------------------------------------------------------
c compute nadapt, ladapt, hadapt 
c---------------------------------------------------------  	  

	  hadapt(ii) = hlv(levmax) 
	  nadapt(ii) = max(0,levmax-lmin(ii)+1) !+1 (min 2 adaptive layer)
	  ladapt(ii) = levmax

c---------------------------------------------------------
c compute hdl: different strategy tested
c--------------------------------------------------------- 

    	  if (nadapt(ii).gt.0) then	  
c	  den = (nsigma(ii)-1.)+r		!freezed
          den = hadapt(ii)-hlv(lmin(ii)-1) !zstar
	  if (ladapt(ii).eq.lmax) den = htot(ii)-hlv(lmin(ii)-1)
          do l=lmin(ii),ladapt(ii)
c	    hdl(l,ii) = - 1. / den		!freezed		
c           hdl(l,ii) = - 1. / nsigma(ii)       !constant
	    htop = hlv(l-1)
	    hbot = hlv(l)
	    if (l.eq.lmax) hbot = htot(ii)
            hdl(l,ii) = (htop-hbot)/den   	!zstar
	  end do
c         hdl(lmin(ii),ii) = - r / den          !freezed
	  check = 0.
	  do l=lmin(ii),ladapt(ii)
            check = check + hdl(l,ii)
          end do		    
	  if (abs(check+1.).gt.1e-5) then
	    write(6,*) 'error computing layer thickness'
	    write(6,*) 'you are using z-star levels'
	    write(6,*) 'but the weights does not sum to -1 ', check
	    write(6,*) 'lmin,ladapt,lmax ', lmin(ii),ladapt(ii),lmax
            write(6,*) 'hadapt,htot ', hadapt(ii),htot(ii)	    
	    stop 'error stop in compute_zadaptive_info'	    
	  end if
	  end if

c---------------------------------------------------------
c compute element info: number of adaptive layers in ele
c---------------------------------------------------------

	  !flag adaptive layers in element with non-conformal edge:
	  !free-surface must span all layers greater then lmin
	  if (ladapt(ii).gt.lmine) then
	    ladapt(4) = max(ladapt(ii),ladapt(4))
	    nadapt(4) = ladapt(4)-lmine+1
	    hadapt(4) = max(hadapt(ii),hadapt(4))	  
	  end if
	end do

	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

        subroutine init_zadaptation

        use levels, only : nlvdi
        use basin, only : nkn,nel		
	use zadapt

	implicit none		

        allocate(nadapt_com(4,nel))
        allocate(hadapt_com(4,nel))	

        nadapt_com  = 0
	hadapt_com  = 0.	
	rzmov_com   = 0

	end

c******************************************************************
c	program zadapt_main
c	call zadapt_test
c	end
c******************************************************************

