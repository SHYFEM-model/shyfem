!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   InitializeSet
!
! FILE
!   InitializeSet.f90
!
! DESCRIPTION
!    Initialize  of all variable in structure  which describe
!     a vertical gradient.
!
! AUTHORS
!    Piet Ruardij, NIOZ
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
      integer function InitializeSet(NUTR,option,input)
        USE global_mem, ONLY:RLEN,ALLOC,error_msg_prn
        use mem, ONLY:NO_BOXES_XY
        USE bennut_variables, ONLY:nutr_seq,nutr_control,ns,sets,fflag
        USE constants, ONLY:LAYERS,NUMBER_OF_PROFILES

        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification: dummy not used!
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification

        integer ::i, NUTo

        if (NUTR == 0 ) then
          if ( nutr_control <=0) then
            allocate(fflag(NUMBER_OF_PROFILES * NO_BOXES_XY), stat=i)
            if (i /= 0) call error_msg_prn(ALLOC,"AllocateMem", "fflag")
            allocate(sets(NUMBER_OF_PROFILES * NO_BOXES_XY), stat=i)
            if (i /= 0) call error_msg_prn(ALLOC,"AllocateMem", "Sets")
            nutr_control=1
            NUTo=nutr_control
          elseif ( nutr_control > NUMBER_OF_PROFILES * NO_BOXES_XY ) then
            stop "InitializeSet:number used profiles>NUMBER_OF_PROFILES"
          else
            nutr_control=nutr_control+1
            NUTo=nutr_control
          endif
        else
          NUTo=NUTR
        endif
        if (option > 0) then
          ns => sets(NUTo)
          nutr_seq=NUTo
          ns%equa=option
          ns%nn=input
          ns%status=LAYERS
          ns%b=-1.D+00
          ns%diff=-1.D+00
          ns%poro=-1.D+00
          ns%ads=-1.D+00
          ns%lst=0
          ns%b(1)=0.0D+00
        endif
        InitializeSet=NUTo
        return
      end

      integer function CopySet(NUTR1,NUTR2)
      USE global_mem, ONLY:RLEN
      USE bennut_variables,ONLY: nutr_control,sets
      USE bennut_interface,ONLY: InitializeSet

      IMPLICIT  NONE
      integer,intent(IN) ::NUTR1 ! Specification
      integer,intent(IN) ::NUTR2 ! Specification

      INTEGER :: idummy =0
      INTEGER :: NUTR2o =0

      ! get  only a new unused sequence number of the list ( if NUTR2 ==0 )
      NUTR2o=InitializeSet(NUTR2,0,idummy)
      sets(NUTR2o)=sets(NUTR1)
      CopySet=NUTR2o
      return
      END

