c
c $Id: ousvdiff.f,v 1.2 2009-04-07 10:43:57 georg Exp $
c
c interpolation of velocities onto nodes
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 04.03.2005	ggu	computes 3D velocities
c 03.11.2010	ggu	prepared to use for water level differences
c 20.12.2014	ggu	bug fix: nlv was not set, new calling sequence
c
c***************************************************************

	program ousdiff

c reads ous files and computes difference

	use mod_depth
	use levels
	use basin

	implicit none

        include 'param.h'

	include 'simul.h'



        real utln1v(nlvdim,neldim)
        real vtln1v(nlvdim,neldim)
        real utln2v(nlvdim,neldim)
        real vtln2v(nlvdim,neldim)

	real weight(nlvdim,nkndim)

	real zn1v(nkndim)
	real zen1v(3,neldim)
	real zn2v(nkndim)
	real zen2v(3,neldim)

	character*80 file1,file2
	integer n,nx,ny
	integer nfreq
	integer ii,l,lmax,nbout,nb1,nb2
	integer ifileo

	logical bnew
        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i,k
	integer it1,it2
        integer ierr,nread,ndry
        integer nknous,nelous,nlvous1,nlvous2
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real cmin,cmax

	integer iapini,ideffi

c--------------------------------------------------------------
c information and copyright
c--------------------------------------------------------------

        call shyfem_copyright('ousdiff - Difference of OUS files')

	write(6,*) 'This program computes the difference of 2 ous files.'
	write(6,*) 'It needs the name of the basin and the 2 ous files.'
	write(6,*) 'The first is subtracted from the second.'

c--------------------------------------------------------------
c open basin
c--------------------------------------------------------------

	nread=0

	if(iapini(1,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

c--------------------------------------------------------------
c open first ous file
c--------------------------------------------------------------

	write(6,*) 'Enter name of first OUS file: '
	read (5,'(a)') file1

        call open_ous_file(file1,'old',nb1)
        call read_ous_header(nb1,nkndim,neldim,nlvdim,ilhv,hlv,hev)
        call ous_get_params(nb1,nknous,nelous,nlvous1)

	if( nkn .ne. nknous ) goto 88
	if( nel .ne. nelous ) goto 88

c--------------------------------------------------------------
c open second ous file
c--------------------------------------------------------------

	write(6,*) 'Enter name of second OUS file: '
	read (5,'(a)') file2

        call open_ous_file(file2,'old',nb2)
        call read_ous_header(nb2,nkndim,neldim,nlvdim,ilhv,hlv,hev)
        call ous_get_params(nb2,nknous,nelous,nlvous2)

	if( nkn .ne. nknous ) goto 88
	if( nel .ne. nelous ) goto 88

c--------------------------------------------------------------
c check nlv
c--------------------------------------------------------------

	if( nlvous1 .ne. nlvous2 ) goto 87
	nlv = nlvous1

        call init_sigma_info(nlv,hlv)

c--------------------------------------------------------------
c open output file
c--------------------------------------------------------------

        call open_ous_file('diff.ous','new',nbout)
        call ous_init(nbout,0)
        call ous_clone_params(nb1,nbout)
        call write_ous_header(nbout,ilhv,hlv,hev)

	write(6,*) 'output file diff.ous opened ',ierr

c--------------------------------------------------------------
c loop on data
c--------------------------------------------------------------

	it1 = -1024*1024*1024
	it2 = it1

  300   continue

	bnew = .true.
	do while( bnew .or. it1 .ne. it2 )
	  do while( bnew .or. it1 .lt. it2 )
            call rdous(nb1,it1,nlvdim,ilhv
     +			,zn1v,zen1v,utln1v,vtln1v,ierr)
	    if( ierr .gt. 0 ) goto 96	!error
	    if( ierr .lt. 0 ) goto 100	!EOF
	    bnew = .false.
	  end do
	  do while( bnew .or. it2 .lt. it1 )
            call rdous(nb2,it2,nlvdim,ilhv
     +			,zn2v,zen2v,utln2v,vtln2v,ierr)
	    if( ierr .gt. 0 ) goto 96	!error
	    if( ierr .lt. 0 ) goto 100	!EOF
	    bnew = .false.
	  end do
	end do

	nread=nread+1
	write(6,*) nread,it1,it2

	do k=1,nkn
	  !zn1v(k) = 0.5 * ( zn2v(k) + zn1v(k) )
	  zn1v(k) = zn2v(k) - zn1v(k)
	end do

	do ie=1,nel
	  do ii=1,3
	    !zen1v(ii,ie) = 0.5 * ( zen2v(ii,ie) + zen1v(ii,ie) )
	    zen1v(ii,ie) = zen2v(ii,ie) - zen1v(ii,ie)
	  end do
	  lmax = ilhv(ie)
	  do l=1,lmax
	    utln1v(l,ie) = utln2v(l,ie) - utln1v(l,ie)
	    vtln1v(l,ie) = vtln2v(l,ie) - vtln1v(l,ie)
	  end do
	end do

        call wrous(nbout,it1,nlvdim,ilhv,zn1v,zen1v,utln1v,vtln1v,ierr)
        if(ierr.ne.0.) goto 95

	goto 300

c--------------------------------------------------------------
c end of loop
c--------------------------------------------------------------

  100	continue

	write(6,*)
	write(6,*) nread,' records elaborated'
	write(6,*)

c--------------------------------------------------------------
c end of routine
c--------------------------------------------------------------

	stop
   87	continue
	write(6,*) 'nlv different between simulations: ',nlvous1,nlvous2
	stop 'error stop: nlv parameter'
   88	continue
	write(6,*) 'nkn, nel from basin: ',nkn,nel
	write(6,*) 'nkn, nel from simul: ',nknous,nelous
	stop 'error stop: simul parameters'
   96	continue
	write(6,*) 'error in reading file : ',ierr
	stop 'error stop: Cannot read file OUS'
   95	continue
	stop 'error stop: Cannot write output file OUS'
	end

c******************************************************************

