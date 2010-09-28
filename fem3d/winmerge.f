c
c $Id: wininf.f,v 1.8 2009-05-29 15:40:42 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c 13.04.1999	ggu	special output introduced
c 08.11.2005	ggu	read also with pressure data
c 07.05.2009	ggu	writes also info on pressure
c 08.09.2010	ggu	new program winmerge
c
c*******************************************************************

	program winmerge

c merges 2 wind files 

c if two same time records are present, 
c the one from the second file is used

	implicit none

        include 'param.h'

	character*80 infile
	real wx1(nkndim)
	real wy1(nkndim)
	real p1(nkndim)
	real wx2(nkndim)
	real wy2(nkndim)
	real p2(nkndim)

	integer ndim
	integer it1,nkn1
	integer it2,nkn2
	integer nin1,nin2,nout
	integer itot
	logical bpres1,bpres2

c---------------------------------------------------------
c initialization
c---------------------------------------------------------

	itot=0
	ndim = nkndim
	nin1 = 1
	nin2 = 2
	nout = 3

c---------------------------------------------------------
c get file name and time shift
c---------------------------------------------------------

	write(6,*) 'Enter name of first wind file :'
	read(5,'(a)') infile
	if(infile.eq.' ') stop
	open(nin1,file=infile,form='unformatted',status='old')

	write(6,*) 'Enter name of second wind file :'
	read(5,'(a)') infile
	if(infile.eq.' ') stop
	open(nin2,file=infile,form='unformatted',status='old')

c---------------------------------------------------------
c loop on wind file
c---------------------------------------------------------

	call read_win(nin1,ndim,it1,nkn1,wx1,wy1,p1,bpres1)
	call read_win(nin2,ndim,it2,nkn2,wx2,wy2,p2,bpres2)

	if( nkn1 .ne. nkn2 ) goto 97

	do while(.true.)

	  if( nkn1 .le. 0 .and. nkn2 .le. 0 ) goto 2
	  itot=itot+1

	  if( nkn1 .gt. 0 .and. it1 .lt. it2 .or. nkn2 .eq. 0 ) then
	    write(6,1000) 1,itot,it1,it1,it2,nkn1,nkn2
	    call write_win(nout,it1,nkn1,wx1,wy1,p1,bpres1)
	    call read_win(nin1,ndim,it1,nkn1,wx1,wy1,p1,bpres1)
	  else if( nkn2 .gt. 0 .and. it2. lt. it1 .or. nkn1 .eq. 0 ) then
	    write(6,1000) 2,itot,it2,it1,it2,nkn1,nkn2
	    call write_win(nout,it2,nkn2,wx2,wy2,p2,bpres2)
	    call read_win(nin2,ndim,it2,nkn2,wx2,wy2,p2,bpres2)
	  else	!equal
	    write(6,1000) 3,itot,it2,it1,it2,nkn1,nkn2
	    call write_win(nout,it2,nkn2,wx2,wy2,p2,bpres2)
	    call read_win(nin2,ndim,it2,nkn2,wx2,wy2,p2,bpres2)
	    call read_win(nin1,ndim,it1,nkn1,wx1,wy1,p1,bpres1)
	  end if
	if( itot .gt. 100 ) stop
	end do

    2	continue
	write(6,*) itot,' wind records written to unit ',nout

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	stop
 1000	format(7i10)
   97	continue
	write(6,*) 'nkn1,nkn2: ',nkn1,nkn2
	stop 'error stop: number of nodes must be the same'
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

