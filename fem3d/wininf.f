c
c $Id: wininf.f,v 1.8 2009-05-29 15:40:42 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c 13.04.1999	ggu	special output introduced
c 08.11.2005	ggu	read also with pressure data
c 07.05.2009	ggu	writes also info on pressure
c
c*******************************************************************

	program wininf

c reads wind file for info

	implicit none

        include 'param.h'

	character*80 infile
	real wx(nkndim)
	real wy(nkndim)
	real s(nkndim)
	real p(nkndim)

	real wxmin,wxmax,wymin,wymax
	real smin,smax
	real pmin,pmax
	integer it,nkn,i
	integer itot
	integer node
	logical bpres

c---------------------------------------------------------------
	node = -1		! -1 for no special node output
	node = 1620
	node = 100
c---------------------------------------------------------------

	wxmin=0.
	wxmax=0.
	wymin=0.
	wymax=0.
	smin=0.
	smax=0.
	itot=0

        do i=1,nkndim
          p(i) = 0.
        end do

	write(6,*) 'Enter name of wind file :'
	read(5,'(a)') infile
	if(infile.eq.' ') stop

	open(1,file=infile,form='unformatted',status='old')

	write(6,*) '   time     nkn  wxmin  wxmax  wymin  wymax  speed'
     +                          ,'     pmin     pmax'

	do while(.true.)
	  read(1,end=2) it,nkn
	  bpres = .false.
	  if( nkn .lt. 0 ) then
	    nkn = -nkn
	    bpres = .true.
	  end if
	  if( nkn .gt. nkndim ) goto 99
          if( bpres ) then
	    read(1) (wx(i),wy(i),i=1,nkn),(p(i),i=1,nkn)
          else
	    read(1) (wx(i),wy(i),i=1,nkn)
          end if
	  do i=1,nkn
	    s(i)=wx(i)**2+wy(i)**2
	  end do
	  call mima(wx,nkn,wxmin,wxmax)
	  call mima(wy,nkn,wymin,wymax)
	  call mima(s,nkn,smin,smax)
	  call mima(p,nkn,pmin,pmax)
	  write(6,'(2i8,5f7.2,2f9.1)') it,nkn
     +			,wxmin,wxmax,wymin,wymax,sqrt(smax)
     +			,pmin,pmax
	  if( node .gt. 0 ) call special(it,node,wx,wy)
	  itot=itot+1
	end do

    2	continue

	close(1)

	write(6,*) itot,' wind records read'

	stop
   99	continue
	write(6,*) 'nkn,nkndim: ',nkn,nkndim
	stop 'error stop wininf: nkn or nkndim'
	end

c*********************************************************************

	subroutine c2p(u,v,s,d)

c converts cartesian to polar coordinates

	implicit none

	real u,v
	real s,d

	real rad,a

c----------------------------------------------------------------
c Pi
c----------------------------------------------------------------

	rad = 45. / atan (1.)

c----------------------------------------------------------------
c compute speed
c----------------------------------------------------------------

	s = sqrt( u**2 + v**2 )

c----------------------------------------------------------------
c compute mathematical angle
c----------------------------------------------------------------

	if( u .eq. 0. ) then
	  if( v .gt. 0. ) then
	    a = 90.
	  else
	    a = -90.
	  end if
	else
	  a = rad * atan(v/u)
	end if

	if( u .lt. 0. ) a = a + 180.

c----------------------------------------------------------------
c convert to wind directions
c----------------------------------------------------------------

	d = 270 - a

	if( d .gt. 360. ) d = d - 360.
	if( d .lt. 0.   ) d = d + 360.

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*********************************************************************

	subroutine special(it,node,wx,wy)

c special output

	implicit none

	integer it
	integer node
	real wx(1), wy(1)

	real u,v,s,d

	u = wx(node)
	v = wy(node)

	call c2p(u,v,s,d)

	write(2,*) it,u,v,s,d

	end

c*********************************************************************

