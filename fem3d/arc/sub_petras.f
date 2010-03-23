
c*****************************************************************
c newini.f
c*****************************************************************

	subroutine inicfil_petras_ali(name,var,nvar)

c initializes nodal value variable from file

	implicit none

	include 'param.h'

	character*(*) name          !name of variable
	real var(nlvdim,nkndim,1)   !variable to set
	real varval
        integer nvar
        integer ftype               !type of  1- homogeneous initial cond.
                                    !2-heterogeneous initial conditions. Added by Petras 12-12-2004 

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	real hlv(1)       !hlv(i)absolute depth of bottom of layer i
	common /hlv/hlv

	character*80 file
	integer nb,irec
	integer nkk,lmax
	integer l,k
        integer ivars,ivar
	real val
	real rlaux(nlvdim) !absolute depth of bottom of layer i
	                   !nlvdim - total number of vertical levels

	integer ifileo

c-------------------------------------------------------
c get file name
c-------------------------------------------------------

	call getfnm(name,file)

	if( file .eq. ' ' ) return	!nothing to initialize

c-------------------------------------------------------
c open file
c-------------------------------------------------------

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value 
C     files for the nodes are read. Before this modification it was not 
C     exactly clear how to format this file. By this modification
C     the file will be read formatted 
C     
C     nb = ifileo(55,file,'unform','old') 
C      -------> nb = ifileo(55,file,'f','old')
C
	nb = ifileo(55,file,'f','old')
	if( nb .le. 0 ) goto 97

c-------------------------------------------------------
c read first record
c-------------------------------------------------------

      ivar = 0

	irec = 1
	
C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value 
C     files for the nodes are read. Before this modification it was not 
C     exactly clear how to format this file. By this modification
C     the file will be read formatted 
C
C     read(nb,err=90) nkk,lmax,ivars 
C      -----> read(nb, 5000, err=90) nkk,lmax,ivars
C
C Skip 4 lines. Added by Petras 12-12-2004     
	read(nb,*)
	read(nb,*)
	read(nb,*)
	read(nb,*)
C Read control information 
	read(nb, 5010, err=90) nkk,lmax,ivars,ftype
 5010 FORMAT(4I5)

      if( nkk .ne. nkn .or. lmax .gt. nlvdim ) goto 99
      if( ivars .ne. nvar ) goto 96
      if(ftype .ne. 1 .and. ftype .ne. 2) goto 91 

C********************************************************
C FILE TYPE 2 (spatialy heterogeneous initial conditions)
C********************************************************
      if(ftype .eq. 1) goto 1000       
c-------------------------------------------------------
c read second record (only if lmax > 0)
c-------------------------------------------------------             !

	if( lmax .gt. 1 ) then          !changed from 0 to 1 (5.3.2004) LMAX
	  irec = 2
C       MODIFIED BY ALI AND PETRAS
C       CORPI, 16 July 2004
C       This modification is related to change the way how inital value 
C       files for the nodes are read. Before this modification it was not 
C       exactly clear how to format this file. By this modification
C       the file will be read formatted 
C	  
C	  
C       read(nb, err=90) (rlaux(l),l=1,lmax)
C       -----> read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
C
	  read(nb, 5020, err=90) (rlaux(l),l=1,lmax)
 5020   FORMAT(8F10.0)
	  
	  do l=1,lmax
	    if( hlv(l) .ne. rlaux(l) ) goto 98
	  end do
	else
	  lmax = 1
	end if

c-------------------------------------------------------
c read data records 
c-------------------------------------------------------

	irec = 3
      
	do ivar=1,nvar

C     MODIFIED BY ALI AND PETRAS
C     CORPI, 16 July 2004
C     This modification is related to change the way how inital value 
C     files for the nodes are read. Before this modification it was not 
C     exactly clear how to format this file. By this modification
C     the file will be read formatted 
C
C     read(nb, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn) 
C --> read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
C     	
	
	read(nb, 5030, err=90) ((var(l,k,ivar),l=1,lmax),k=1,nkn)
 5030 FORMAT(F5.0)

      end do
      goto 1001
      
C******************************************************
C FILE TYPE 1 (spatialy homogeneous initial conditions)
C******************************************************
 1000 continue
      irec=3
      do ivar=1,nvar
       read(nb, 5030, err=90) varval 
       do l=1,lmax
        do k=1,nkn
         var(l,k,ivar)=varval
        end do
       end do
      end do    

c-------------------------------------------------------
c reading done -> close file
c-------------------------------------------------------
 1001 continue  
	close(nb)

c-------------------------------------------------------
c initialize the other levels if only surface is given
c-------------------------------------------------------

        do ivar=1,nvar
	  do k=1,nkn
	    val = var(lmax,k,ivar)
	    do l=lmax+1,nlvdim
	      var(l,k,ivar) = val
	    end do
	  end do
        end do

c-------------------------------------------------------
c end of routine
c-------------------------------------------------------

	write(6,*) 'Succesfull initialization for ',name,' from file '
	write(6,*) file

	return
   90	continue
	write(6,*) 'read error in record = ',irec,' ivar = ',ivar
	write(6,*) '... reading file',file
	stop 'error stop inicfil'
91	continue
	write(6,*) 'bad file type descriptor: value 1 or 2 is allowed only!' 
	write(6,*) '... reading file',file
	stop 'error stop inicfil'
   96	continue
	write(6,*) 'ivars not compatible with nvar: ',ivars,nvar
	stop 'error stop inicfil'
   97	continue
	write(6,*) 'Cannot open file ',file
	stop 'error stop inicfil'
   98	continue
	write(6,*) 'levels are not the same from init file ',file
	write(6,*) (hlv(l),l=1,lmax)
	write(6,*) (rlaux(l),l=1,lmax)
	stop 'error stop inicfil'
   99	continue
	write(6,*) 'parameters are not the same from init file ',file
	write(6,*) 'nkn, lmax from file  : ',nkk,lmax
	write(6,*) 'nkn, lmax from model : ',nkn,nlvdim
	stop 'error stop inicfil'
	end

c*******************************************************************
c debug.f
c***************************************************************

      subroutine check2Dbio(it,ul,nlv,n,nvar,a,vmin,vmax,textgen,text)

c tests array for nan and strange values for EUTRO variables
c Made from check2Dr by Petras 11 January 2005      

	implicit none
	
	include 'param.h'
	
	integer ipv(nkndim) !external numbers of nodes
	common /ipv/ipv
	integer nlv,n
	integer nvar         !state variable number
	real a(nlvdim,1)     !array of variable values for levels and nodes
	real vmin,vmax       !minimal and maximal allowed values
	character*(*) textgen,text ! '***BIO CHECK' and text indicating the
	                           ! time of checking

	logical debug,bval
	integer inan,iout,i,l,ul,it
	real val

	logical is_r_nan

        bval = vmin .lt. vmax
	debug = .true.
	inan = 0
	iout = 0

	do i=1,n
	  do l=1,nlv
	    val = a(l,i)	    
	    if( is_r_nan(val) ) then
	      inan = inan + 1
	      if( debug ) write(ul,'(I2,G11.3,I5,I2,I10,1X,A16)') nvar,val,
     +                        ipv(i),l,it,text	      
	    else if( bval .and. (val .lt. vmin .or. val .gt. vmax) ) then
	      iout = iout + 1	      
	      if( debug ) write(ul,'(I2,G11.3,I5,I2,I10,1x,A16)') nvar,val,
     +                        ipv(i),l,it,text
     	    end if
	  end do
	end do

	if( inan .gt. 0 .or. iout .gt. 0 ) then
	  write(6,'(2A13,A2,A16,A2,2(A5,I4))') 'CHECK2DBIO: ',
     +	textgen," (",text,") ",
     +  ' NaN=',inan,' OUT=',iout	  
	end if
	
	end
	  
c***************************************************************


