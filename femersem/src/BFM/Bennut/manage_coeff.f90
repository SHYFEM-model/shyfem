!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
!
! FUNCTION
!   store_coeff.f90
!   read_coeff.f90
!   select_coeff.f90
!
! FILE
!   manage_coeff.f90
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
      integer FUNCTION store_coeff(file,box)
        USE global_mem, ONLY:RLEN
        USE bennut_interface,ONLY:io_coeff
        USE bennut_constants
        IMPLICIT  NONE
        character file*(*) ! Specification
        integer,intent(IN) ::box ! Specification

        integer ::start
        integer ::i
        integer ::l
        integer ::k
        integer ::first
        save start,first
        REAL(RLEN) ::idummy

        character hulp*(26)

        data start/1/,first/0/

        if (first.eq.0) then
          first=box
          k=1
        elseif( first.eq.box) then
          k=1
          start=start+1
        else
          k=k+1
        endif
        if (start.eq.1) then
          i=index(file,'.')
          l=len_trim(file)
          if (i.gt.0) l=i-1
          write(hulp(l:),'(I5.5)')box
          hulp(1:l)=file(1:l)
          hulp(l+5:)=file(l+1:)
          open(unit=200+k,file=hulp,access='sequential', &
                                               form='unformatted')
        endif

        if (start.eq.1) then
           i=box
        else
           i=-box
        endif
        idummy= io_coeff(WRITE_COEFF,200+k,i)

        store_coeff=start

        return
        end

!
      integer FUNCTION read_coeff(file)
        USE global_mem, ONLY:RLEN
        USE bennut_interface,ONLY:io_coeff
        character,intent(IN) ::file*(*) ! Specification
        

        integer ::start
        integer ::idummy
        save start

        data start/-1/


        REAL(RLEN) ::i

        if (start.eq.-1) then
          open(unit=93,file=file,access='sequential',form='unformatted')
        endif

        i= io_coeff(READ_COEFF,93,idummy)
        start=start-1

        if (i.gt.0) then
          read_coeff=-start
        else
          read_coeff=-1
        endif

        return
      end

     integer FUNCTION select_coeff(mode,box)
        USE global_mem, ONLY:RLEN
        USE bennut_interface,ONLY:io_coeff
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::box ! Specification

        REAL(RLEN) ::i


        i= io_coeff(SELECT_COEFF,mode,box)

        select_coeff=i

        return
      end


      integer FUNCTION io_coeff(mode,option,input) 
        USE global_mem, ONLY:RLEN
        USE mem, ONLY:NO_BOXES_XY
        USE bennut_type
        USE bennut_variables
        USE bennut_constants
        USE constants, ONLY:NUMBER_OF_PROFILES
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::option ! Specification
        integer,intent(IN) ::input ! Specification

        integer ::i,nmax
        integer ::j


        nmax=NUMBER_OF_PROFILES * NO_BOXES_XY
        IF (mode == WRITE_COEFF) then
          !write.....
          !write head
          if (input.gt.0) then
            if (nflag < 0) nflag=-nflag
            j=0
            do i=1,nmax
              if (fflag(j) == input) j=j+1
            enddo
            write(option) j,input
          elseif(input == 0) then
            stop 'io_coeff:WRITE_COEFF box=0 && with selection'
          endif
          do j=1,nmax
            if (fflag(j) == abs(input)) write(option) sets(j)
          enddo
        elseif (mode == READ_COEFF) then
          !read....
          !read head
          if (nflag == 0) then
            READ(option) nmax,i
            fflag(1:nmax)=i
            nflag=nmax
          endif
          do j=1,nmax
            read(option) sets(j)
          enddo
        elseif (mode == SELECT_COEFF) then
          if (nflag == 0) fflag(1:nmax)=0
          if (option.gt.nmax) then
            stop &
             'io_coeff: sequence number of set to be saved out of range'
          elseif(nflag.gt.0) then
            write(0,*) 'io_coeff:', &
             'sets for saving only be set before the first writing'
            stop
          else
            nflag=nflag-1
            fflag(option)=input
          endif
        ENDIF
        io_coeff=1
      end
