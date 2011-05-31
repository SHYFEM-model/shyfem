!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   CalculateShift
!
! FILE
!   CalculateShift.f90
!
! DESCRIPTION
!   
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
!   
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!

     real(RLEN) function CalculateShift(KMI,layer,from,to)
     USE global_mem, ONLY:RLEN
     USE bennut_interface,ONLY: noutput,CalculateLayer
     USE constants, ONLY: LAYERS,INTEGRAL,MASS
     implicit none 
     integer    :: layer,KMI       ! Specification
     real(RLEN) ::  from,to        ! Specification

      integer ::  n,i
      real(RLEN)    ::  rshift
      real(RLEN)    ::  value

      ! get number of layers (in n )
      ! value and rchift are used as dummy vars.
      call CalculateLayer(KMI,-1,value,n,rshift)
      rshift=to-from
      i=layer
      if ( layer > n ) then
         stop 'CalculateShift: layer>= nr_layers '
      elseif (layer < n) then
         if ( rshift  > 0.0D+00) i=i+1
      endif 

      ! Calculate C_MASS of nutrients which is rshifted                    
      value = noutput(KMI,i, INTEGRAL, MASS, from, to)

      ! Check Concistency : for very small rshifts * delts it may happen
      ! that the answer has a other sighn than the rshift                
      if (value * rshift < 0.0) value = 0.0

       CalculateShift=value
       return

       end


