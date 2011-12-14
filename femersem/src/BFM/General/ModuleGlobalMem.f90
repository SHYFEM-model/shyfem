!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleGlobalMem
!
! DESCRIPTION
!    !   Definiton of the runtime error messages.
!   This module contains global settings:
!   -general constants for controlling prescision,
!   -parameters defining fle streams and error message numbers
!   -the subroutine for printing the message
!   and aborting the simulation

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE global_mem
!

!  
!
! !AUTHORS
!   P. Ruardij (NIOZ)/ M. Vichi (INGV) 
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Global Constants
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! REALS are double precision
  integer,    parameter ::RLEN=kind(1D0)
  integer,    parameter ::NMLUNIT=310
  integer,    parameter ::LOGUNIT=0
  real(RLEN), parameter ::ZERO=0.0D+00
  real(RLEN), parameter ::PI=3.14159265359D+00
  real(RLEN), parameter ::BASETEMP= 10.0D+00
  real(RLEN), parameter ::ZERO_KELVIN=-273.3D+00
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Next parameters are defined to define differnt types of State variables. &
  ! There are four different types defined:
  integer,    parameter ::SINKSOURCE=-1
  integer,    parameter ::NOTRANSPORT=0
  integer,    parameter ::HORTRANSPORT=10
  integer,    parameter ::ALLTRANSPORT=20
  real(RLEN), parameter :: DONE=1.D+00
  ! Error codes:
  integer,    parameter ::ALLOC=10
  integer,    parameter ::NML_OPEN=11
  integer,    parameter ::NML_READ=12
  integer,    parameter ::DIM_MISMATCH=13
  contains

  subroutine error_msg_prn(code,infile,what)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer :: code
  character(LEN=*) :: what,infile
write(LOGUNIT,*) "*********** RUN TIME ERROR BEGIN ***********"
      select case (code)
        case (ALLOC)
          write(LOGUNIT,*) "Unable to allocate ",trim(what)," in ", &
                                trim(infile)
        case (NML_OPEN)
          write(LOGUNIT,*) "Unable to open ",trim(what)," in ",  &
                                trim(infile)
        case (NML_READ)
          write(LOGUNIT,*) "Namelist mismatch in ",trim(what),  &
                         " opened by ",trim(infile)
        case (DIM_MISMATCH)
          write(LOGUNIT,*) "Dimension mismatch while reading ",  &
                trim(what)," in ",trim(infile)
    end select
    write(LOGUNIT,*) "***********  RUN TIME ERROR END  ***********"
    stop "SIMULATION ABORTED (see MEM logfile)"
  end subroutine
  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
