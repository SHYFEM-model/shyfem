c
c $Id: plotsim.f,v 1.53 2010-03-11 15:33:09 georg Exp $
c
c revision log :
c
c 12.02.1999	ggu	adapted to auto mode
c 30.10.2003    ggu     added plowind()
c 22.03.2004	ggu	lagrang model output added (mode 12)
c 04.10.2004	ggu	temp and salt from nos file
c 02.12.2004    ggu     copyright statement added
c 02.03.2005    ggu     introduced flag for scalar plot
c 18.10.2005    ggu     data structures introduced to call set_geom
c 14.03.2007    ggu     wave plotting
c 16.04.2008    ggu     new Makefile structure
c 09.12.2008    ggu     changes from Malta integrated (annotation)
c 26.01.2009	ggu	new makedepend
c 06.04.2009	ggu	new param.h structure
c 20.04.2009	ggu	level.h eliminated
c 09.10.2009	ggu	plotting of atmos. pressure
c 13.10.2009	ggu	new arrays for velocity section plot
c 23.02.2010	ggu	new call to set_default_color_table()
c 17.12.2010	ggu	substituted hv with hkv
c 31.03.2011	ggu	new arrays fvlv, arfvlv for scalar plotting
c 31.08.2011	ggu	new copyright, eos plotting
c 14.11.2011	ggu	new array hl for depth structure
c
c*************************************************************

	program plotsim

c plots simulation

	implicit none

c parameters

	include 'param.h'
	include 'basin.h'
	include 'evmain.h'

	integer matdim
	parameter (matdim=nkndim*ngrdim)

c description
	character*80 descrr
	common /descrr/ descrr
c FEM parameters
	real eps,eps2,pi,flag,high
	common /mkonst/ eps,eps2,pi,flag,high
c OUT
	real xv(3,nkndim)
	common /xv/xv
	real zenv(3,neldim)
	common /zenv/zenv
	real usnv(neldim), vsnv(neldim)
	common /usnv/usnv, /vsnv/vsnv
	real wsnv(nkndim)
	common /wsnv/wsnv
c 3d
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real hlv(nlvdim), hev(neldim)
        common /hlv/hlv, /hev/hev
        integer ilhv(neldim)
        common /ilhv/ilhv
        integer ilhkv(nkndim)
        common /ilhkv/ilhkv
        real znv(nkndim)
        common /znv/znv
        real utlnv(nlvdim,neldim)
        common /utlnv/utlnv
        real vtlnv(nlvdim,neldim)
        common /vtlnv/vtlnv
        real ulnv(nlvdim,neldim)
        common /ulnv/ulnv
        real vlnv(nlvdim,neldim)
        common /vlnv/vlnv
        real het3v(nlvdim,neldim)
        common /het3v/het3v

        real uprv(nlvdim,nkndim)
        common /uprv/uprv
        real vprv(nlvdim,nkndim)
        common /vprv/vprv
        real wprv(nlvdim,nkndim)
        common /wprv/wprv

        real p3(nlvdim,neldim)		!is good for nodes and elements
        common /p3/p3

        real fvlv(nlvdim,nkndim)
        common /fvlv/fvlv
        real arfvlv(nkndim)
        common /arfvlv/arfvlv

c boundary etc.
	integer kantv(2,nkndim)
	common /kantv/kantv
        integer ilinkv(nkndim+1),lenkv(nlkdim),linkv(nlkdim)
        common /ilinkv/ilinkv, /lenkv/lenkv, /linkv/linkv
        integer ieltv(3,neldim)
        common /ieltv/ieltv
	real hkv(nkndim)
	real hetv(neldim)
	real parray(neldim)		!is good for nodes and elements
	common /hkv/hkv
	common /hetv/hetv
	common /parray/parray
c auxiliary
	real v1v(nkndim)
	common /v1v/v1v
	real v2v(nkndim)
	common /v2v/v2v
	real v3v(nkndim)
	common /v3v/v3v
	real vev(neldim)
	common /vev/vev
	real amat(matdim)
	common /amat/amat
	real uvnv(neldim), vvnv(neldim)
	common /uvnv/uvnv, /vvnv/vvnv
	real uv(nkndim), vv(nkndim)
	common /uv/uv, /vv/vv
	real pres(nkndim)
	common /pres/pres
	real hl(nlvdim)
	common /hl/hl
	logical bwater(neldim)
	common /bwater/bwater
	logical bkwater(nkndim)
	common /bkwater/bkwater
c vertical velocity
        real wauxv(0:nlvdim,nkndim)
        common /wauxv/wauxv
        real wlnv(0:nlvdim,nkndim)
        common /wlnv/wlnv
c local
	integer mode
	integer ie,ii,k,l,i
	integer icolor,isphe
	integer iapini
	real sflag
	real getpar

c----------------------------------------------
c copyright
c----------------------------------------------

	call shyfem_copyright('plotsim - plotting maps on FE grids')

c----------------------------------------------
c initialize parameters
c----------------------------------------------

	eps=1.e-5
	eps2=1.e-6
	pi=3.141592653
	flag=-9988765.0
	high=1.e+35
	grav=9.81

	nlvdi = nlvdim
	nlv = nlvdim

	call set_flag(sflag)

	do i=1,neldim
	  do l=1,nlvdim
	    p3(l,i) = sflag
	  end do
	end do

c----------------------------------------------
c read basin
c----------------------------------------------

	if(iapini(7,nkndim,neldim,matdim).eq.0) then
		stop 'error stop : iapini'
	end if

c----------------------------------------------
c set up elements
c----------------------------------------------

	call setdim('nlkdim',nlkdim)
	isphe = nint(getpar('isphe'))
	call set_coords_ev(isphe)
	call set_ev
	call set_geom

c----------------------------------------------
c make depth on nodes and elements
c----------------------------------------------

	call mkhkv(hkv,v1v,nkn,nel)
	call mkhev(hev,nel)
	call mkareafvl			!area of finite volumes

c----------------------------------------------
c interactive set up
c----------------------------------------------

	call colsetup
	icolor = nint(getpar('icolor'))
	call set_color_table( icolor )
	call set_default_color_table( icolor )


	call asklev		!ask for 3d level

	call ichoice(mode)

c----------------------------------------------
c open plot
c----------------------------------------------

	call qopen

c----------------------------------------------
c do plotting
c----------------------------------------------

	if( mode .eq. 1 )  call plobas
	if( mode .eq. 2 )  call plosim(.true.)
	if( mode .eq. 3 )  call plosim(.false.)
	if( mode .eq. 4 )  call plozet
c	if( mode .eq. 4 )  call plobar
	if( mode .eq. 5 )  call plonos('.con',10)
	if( mode .eq. 6 )  call plonos('.nos',12)
	if( mode .eq. 7 )  call plonos('.nos',11)
	if( mode .eq. 8 )  call plonos('.rms',18)
	if( mode .eq. 9 )  call plonos('.oxy',15)
	if( mode .eq. 10 ) call plonos('.nos',0)
	if( mode .eq. 11 ) call plowind
	if( mode .eq. 12 ) call plolagr
	if( mode .eq. 13 ) call plowave
	if( mode .eq. 14 ) call plopres
	if( mode .eq. 15 ) call ploeos('.eos',0)

c----------------------------------------------
c close plot
c----------------------------------------------

	call qclose

c----------------------------------------------
c end of routine
c----------------------------------------------

	end

c*****************************************************************

	subroutine ichoice(mode)

	implicit none

	integer mode

	integer iwhat,iauto
	integer ideflt
	real getpar

	iwhat = nint(getpar('iwhat'))
	iauto = nint(getpar('iauto'))

	write(6,*)
	write(6,*) ' basin ...................  1'
	write(6,*) ' velocity ................  2'
	write(6,*) ' transport ...............  3'
	write(6,*) ' water level .............  4'
	write(6,*) ' concentration ...........  5'
	write(6,*) ' temperature .............  6'
	write(6,*) ' salinity ................  7'
	write(6,*) ' rms .....................  8'
	write(6,*) ' oxygen ..................  9'
	write(6,*) ' generic concentrations .. 10'
	write(6,*) ' wind array .............. 11'
	write(6,*) ' lagrangian .............. 12'
	write(6,*) ' wave ..... .............. 13'
	write(6,*) ' atmos. pressure ......... 14'
	write(6,*) ' generic element values .. 15'
	write(6,*)

        if( iauto .eq. 0 .or. iwhat .eq. 0 ) then
          iwhat = ideflt(iwhat,'Enter choice : ')
        else
          write(6,*) 'Plotting : ',iwhat
	  write(6,*)
        end if

	mode = iwhat

	end

c*****************************************************************

