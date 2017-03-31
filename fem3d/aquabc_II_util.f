! AQUABC utilities independently of version
! Contains:
! subroutine set_3d_r_array
! subroutine set_2d_r_array
! subroutine set_1d_r_array
! subroutine set_3d_d_array
! subroutine set_2d_d_array
! subroutine set_1d_d_array
! subroutine dump_aquabc
! subroutine fluxes_aquabc
!
!

c********************************************************************
c version for reals
c********************************************************************

	subroutine set_3d_r_array(n1,n2,n3,ra,rval)

	implicit none

	integer n1,n2,n3
	real ra(n1,n2,n3)
	real rval

	integer i1,i2,i3

	do i3=1,n3
	  do i2=1,n2
	    do i1=1,n1
	      ra(i1,i2,i3) = rval
	    end do
	  end do
	end do

	end

c********************************************************************

	subroutine set_2d_r_array(n1,n2,ra,rval)

	implicit none

	integer n1,n2
	real ra(n1,n2)
	real rval

	integer i1,i2

	do i2=1,n2
	  do i1=1,n1
	    ra(i1,i2) = rval
	  end do
	end do

	end

c********************************************************************

	subroutine set_1d_r_array(n,ra,rval)

	implicit none

	integer n
	real ra(n)
	real rval

	integer i

	do i=1,n
	  ra(i) = rval
	end do

	end

c********************************************************************
c version for double precision
c********************************************************************

	subroutine set_3d_d_array(n1,n2,n3,da,dval)

	implicit none

	integer n1,n2,n3
	double precision da(n1,n2,n3)
	double precision dval

	integer i1,i2,i3

	do i3=1,n3
	  do i2=1,n2
	    do i1=1,n1
	      da(i1,i2,i3) = dval
	    end do
	  end do
	end do

	end

c********************************************************************

	subroutine set_2d_d_array(n1,n2,da,dval)

	implicit none

	integer n1,n2
	double precision da(n1,n2)
	double precision dval

	integer i1,i2

	do i2=1,n2
	  do i1=1,n1
	    da(i1,i2) = dval
	  end do
	end do

	end

c********************************************************************

	subroutine set_1d_d_array(n,da,dval)

	implicit none

	integer n
	double precision da(n)
	double precision dval

	integer i

	do i=1,n
	  da(i) = dval
	end do

	end

c********************************************************************
c********************************************************************
c   Dumps state variables for repeated runs   
c
       subroutine dump_aquabc(file_name,state, 
     +                 nkndim,nkn,nlvdim,nvar)
     
      implicit none
      character*(*) file_name
      integer nkndim,nkn,nlvdim,nvar
      integer un
      real state(nlvdim,nkndim,nvar)
      
      integer NUM_COLS, FILE_TYPE
      character * 5  NUM_STRING
      character * 30 FORMAT_STRING
      integer ifileo
      integer i,j,k
      
      NUM_COLS = nlvdim
      write(NUM_STRING,100) NUM_COLS
      FORMAT_STRING = '(' // NUM_STRING // 'F20.8)'
      
      un = ifileo(55,file_name,'form','unknown')
      
      write(un,*) '***************************************'
      write(un,*) '* DUMP FILE CREATED BY SHYFEM-AQUABC  *'
      write(un,*) '* !!! PLEASE DO NOT EDIT MANUALLY !!! *'
      write(un,*) '***************************************'
      
      FILE_TYPE = 4 
      
      write(un,101) nkn, NUM_COLS, nvar, FILE_TYPE
       
           
      do i = 1, nvar
          do j = 1, nkn
              
              write(un,FORMAT_STRING) (state(k,j,i),k=1,NUM_COLS)
          end do
      end do
      
      
  100 format(i5)    
  101 format(4i5)    
      close(un)
      
      return
      end



c********************************************************************
c********************************************************************      
        subroutine fluxes_aquabc(it,iconz,conzv)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Administers writing of flux data for  ecological variables
c 
c 
c 
c 
c 
c 
c 
c Inputs:
c  it     - Time(seconds)
c  iconz  - Actual number of variables that fluxes should be calculated
c  conzv  - Array of variable values that fluxes should be calculated
c 
c Other important variables:
c ncsdim -    Maximum number of variables that fluxes should be calculated, defined in param.h
c csc    -    New extension for fluxes file
c ivar_base - Base of variable numbering in flux file (harcoded 200 in this subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none

        include 'param.h'

        integer it

        integer nscdim
        parameter(nscdim=20)			!maximum number of sections

        integer nsect,kfluxm,kflux(1)
        common /kfluxc/ nsect,kfluxm,kflux
        integer iflux(3,1)
        common /iflux/iflux
        save /kfluxc/,/iflux/	!ggu

        real conzv(nlvdim,nkndim,ncsdim)	!multiple concentrations
!        common /conzv/conzv

        integer itend
        integer j,i,k,l,lmax,nlmax,ivar,nvers,ivar_base
        integer iconz
        real az,azpar,rr,dt

        real fluxes(0:nlvdim,3,nscdim)

        integer nlayers(nscdim)			!number of layers in section
        real trc(ncsdim)			!counter
        real cflux(0:nlvdim,3,nscdim,ncsdim)	!accumulator

        save nlayers,trc,cflux

        integer idtflx,itflx,itmflx,nbflx
        save idtflx,itflx,itmflx,nbflx

        data nbflx /0/

       integer ifemop
       real getpar

c-----------------------------------------------------------------
c start of code
c-----------------------------------------------------------------


        if( nbflx .eq. -1 ) return

c-----------------------------------------------------------------
c initialization
c-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

                idtflx = nint(getpar('idtflx'))
                itmflx = nint(getpar('itmflx'))
                itend = nint(getpar('itend'))

!             iconz = nint(getpar('iconz')) !computing concentrations?

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( idtflx .le. 0 ) nbflx = -1
                if( itmflx .gt. itend ) nbflx = -1
                if( iconz .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

                if( nsect .gt. nscdim ) then
                  stop 'error stop fluxes_template: dimension nscdim'
                end if

                itflx = itmflx + idtflx
                itmflx = itmflx + 1      !start from next time step

                call get_nlayers(kfluxm,kflux,nlayers,nlmax)

                do k=1,iconz
                 call fluxes_init(nlvdim,nsect,nlayers
     +				,trc(k),cflux(0,1,1,k))
                end do

                nbflx=ifemop('.csc','unform','new')
                if(nbflx.le.0) then
                 stop 'error stop wrflxa : Cannot open csc file'
                end if

                nvers = 5
                call wfflx      (nbflx,nvers
     +                          ,nsect,kfluxm,idtflx,nlmax
     +                          ,kflux
     +                          ,nlayers
     +                          )

        end if
                                                  
c-----------------------------------------------------------------
c normal call
c-----------------------------------------------------------------

        if( it .lt. itmflx ) return

!        iconz = nint(getpar('iconz'))

	call get_timestep(dt)
        call getaz(azpar)
        az = azpar
        ivar_base = 200		!base of variable numbering

c	-------------------------------------------------------
c	accumulate results
c	-------------------------------------------------------



        do k=1,iconz
         ivar = ivar_base + k
         call flxscs(kfluxm,kflux,iflux,az,fluxes,ivar,conzv(1,1,k))
         call fluxes_accum(nlvdim,nsect,nlayers
     +			,dt,trc(k),cflux(0,1,1,k),fluxes)
        end do

        
!      print *,'llllllllllllllllllllllllllllllllllllllllllllll'
!      print *, 'iconz= ',iconz
!      stop 
c	-------------------------------------------------------
c	time for output?
c	-------------------------------------------------------

        if( it .lt. itflx ) return
        itflx=itflx+idtflx

c	-------------------------------------------------------
c	average and write results
c	-------------------------------------------------------

        do k=1,iconz
         ivar = ivar_base + k
         call fluxes_aver(nlvdim,nsect,nlayers
     +			,trc(k),cflux(0,1,1,k),fluxes)
         call wrflx(nbflx,it,nlvdim,nsect,ivar,nlayers,fluxes)
        end do

c	-------------------------------------------------------
c	reset variables
c	-------------------------------------------------------

        do k=1,iconz
         call fluxes_init(nlvdim,nsect,nlayers
     +			,trc(k),cflux(0,1,1,k))
        end do

c-----------------------------------------------------------------
c end of routine fluxes_aquabc
c-----------------------------------------------------------------

        end

c******************************************************************
c******************************************************************
