c
c $Id: wininf.f,v 1.8 2009-05-29 15:40:42 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c 13.04.1999	ggu	special output introduced
c 08.11.2005	ggu	read also with pressure data
c 07.05.2009	ggu	writes also info on pressure
c 08.09.2010	ggu	new program winshift
c
c*******************************************************************

	program winshift

c reads wind file and shifts time 

	implicit none

        include 'param.h'

	character*80 infile
	real wx(nkndim)
	real wy(nkndim)
	real p(nkndim)

	integer it,idt,nkn,ndim
	integer nin,nout
	integer itot
	logical bpres

c---------------------------------------------------------
c initialization
c---------------------------------------------------------

	itot=0
	ndim = nkndim
	nin = 1
	nout = 2

c---------------------------------------------------------
c get file name and time shift
c---------------------------------------------------------

	write(6,*) 'Enter name of wind file :'
	read(5,'(a)') infile
	if(infile.eq.' ') stop

	write(6,*) 'Enter time shift :'
	read(5,*) idt
	if(idt.eq.0) stop

	open(nin,file=infile,form='unformatted',status='old')

c---------------------------------------------------------
c loop on wind file
c---------------------------------------------------------

	do while(.true.)
	  call read_win(nin,ndim,it,nkn,wx,wy,p,bpres)
	  if( nkn .le. 0 ) goto 2
	  it = it + idt
	  call write_win(nout,it,nkn,wx,wy,p,bpres)
	  itot=itot+1
	  write(6,*) itot,it,nkn,bpres
	end do

    2	continue
	write(6,*) itot,' wind records written to unit ',nout

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c*********************************************************************

	subroutine write_win(nout,it,nkn,wx,wy,p,bpres)

	implicit none

	logical bpres
	integer nout,it,nkn
	real wx(1), wy(1), p(1)

	integer i

        if( bpres ) then
	  write(nout) it,-nkn
	  write(nout) (wx(i),wy(i),i=1,nkn),(p(i),i=1,nkn)
        else
	  write(nout) it,nkn
	  write(nout) (wx(i),wy(i),i=1,nkn)
        end if

	end

c*********************************************************************

	subroutine read_win(nin,ndim,it,nkn,wx,wy,p,bpres)

	implicit none

	logical bpres
	integer nin,ndim,it,nkn
	real wx(1), wy(1), p(1)

	integer i

	nkn = 0
	read(nin,end=2) it,nkn

	if( nkn .lt. 0 ) then
	  nkn = -nkn
	  bpres = .true.
	end if
	if( nkn .gt. ndim ) goto 99

        if( bpres ) then
	  read(nin) (wx(i),wy(i),i=1,nkn),(p(i),i=1,nkn)
        else
	  read(nin) (wx(i),wy(i),i=1,nkn)
        end if

    2	continue
	return
   99	continue
	stop 'error stop read_win: ndim'
	end

c*********************************************************************

