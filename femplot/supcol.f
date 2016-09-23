c
c $Id: supcol.f,v 1.12 2010-02-26 15:29:19 georg Exp $
c
c color handling routines
c
c revision log :
c
c 26.10.2001    ggu     in colrd allow only one array to be given
c 06.06.2008    ggu     new colortable introduced
c 19.11.2008    ggu     new routine coltst()
c 09.01.2009    ggu     some re-formatting
c 09.01.2009    ggu     deleted setcol(), new make_single_isolines()
c 26.01.2009    ggu     some more color tables
c 23.02.2010    ggu     new set_default_color_table(), renamed qchsc()
c 23.03.2010    ggu     make it easier to change between color tables
c
c************************************************************

	subroutine colini

c initializes color

	use color

	implicit none

	logical, save :: binit = .false.

	if( binit ) return

	binit = .true.

	write(6,*) 'initializing color with colini...'

	isopar = isodim
	isoanz = 0
	nisord = 0
	ncolrd = 0
	icauto = -1
	fnull = -999.0
	ciso(1) = -1.

	end

c************************************************************

	subroutine colrd

c reads color block from str file

	use color

	implicit none

	character*80 name,text
	logical bdebug
	integer mode
	integer i
	integer ndim
	real value
	double precision dvalue

	integer nrdpar
	real getpar

	call colini

	ndim = isodim
	bdebug = .true.
	bdebug = .false.
	ncolrd = 0
	nisord = 0

        mode = 1
        do while( mode .ne. 0 )
            mode = nrdpar('color',name,dvalue,text)
	    value = dvalue
	    !write(6,*) mode,name,value,text
            if( mode .eq. 2 ) then
	      if( name .ne. 'color' .and. name .ne. 'isoval' ) then
		write(6,*) 'Variable is no vector : ',name
		stop 'error stop : colrd'
	      end if
	    end if
            if( name .eq. 'color' ) then
                ncolrd = ncolrd + 1
		if( ncolrd .gt. ndim+1 ) stop 'error stop colrd: dimension'
                ciso(ncolrd) = value
            else if( name .eq. 'isoval' ) then
                nisord = nisord + 1
		if( nisord .gt. ndim ) stop 'error stop colrd: dimension'
                fiso(nisord) = value
	    end if
        end do

	if( ncolrd .gt. 0 .and. nisord .gt. 0 ) then
	  if( nisord + 1 .ne. ncolrd ) then
		write(6,*) 'ncolrd,nisord: ',ncolrd,nisord
		write(6,*) 'There must be one more color than isovalue'
		stop 'error stop : colrd'
	  end if
	end if

	isoanz = nisord
c	call putpar('dval',0.)	!if section is given -> no automatic

	if( bdebug ) then
	  write(6,*) 'color debug : ',nisord,ncolrd
	  write(6,*) (fiso(i),i=1,nisord)
	  write(6,*) (ciso(i),i=1,ncolrd)
	end if

	end

c************************************************************

	subroutine coltst

c debug write of color

	use color

	implicit none

	integer i

	write(6,*) 'coltst: debug write of color common block...'

	write(6,*) 'isodim : ',isodim
	write(6,*) 'isopar : ',isopar
	write(6,*) 'isoanz : ',isoanz
	write(6,*) 'nisord : ',nisord
	write(6,*) 'ncolrd : ',ncolrd
	write(6,*) 'icauto : ',icauto
	write(6,*) 'fnull : ',fnull

	if( nisord .lt. 0 .or. nisord .gt. 1000 ) goto 99
	if( ncolrd .lt. 0 .or. ncolrd .gt. 1000 ) goto 99

	write(6,*) 'fiso : ',(fiso(i),i=1,nisord)
	write(6,*) 'ciso : ',(ciso(i),i=1,ncolrd)

	return
   99	continue
	write(6,*) 'unrealistic values for nisord or ncolrd'
	stop 'error stop coltst: nisord/ncolrd'
	end

c************************************************************
c************************************************************
c************************************************************
c accessor routines
c************************************************************
c************************************************************
c************************************************************

	subroutine colcopy(nval,fisov,col)

c sets colors and isolevels

	use color

	implicit none

	integer nval
	real fisov(1), col(1)
	integer i

	call colini

	if( nval .gt. isopar .or. nval .lt. 0 ) then
		write(6,*) 'nval,isopar: ',nval,isopar
		stop 'error stop colcopy: error in parameters'
	end if

	isoanz = nval

	do i=1,nval
	  fiso(i) = fisov(i)
	  ciso(i) = col(i)
	end do

	ciso(nval+1) = col(nval+1)

	end
 
c************************************************************

	subroutine colnul(fnul)

c sets null value

	use color

	implicit none

	real fnul

	call colini

	fnull = fnul

	end
 
c************************************************************

	subroutine colinfo(numiso,fnul)

c returns info on color table

	use color

	implicit none

	integer numiso
	real fnul

	call colini

	numiso = isoanz
	fnul = fnull

	end
 
c************************************************************

	subroutine colminmax(vmin,vmax)

c returns minimum and maximum value in color table

	use color

	implicit none

	real vmin,vmax

	call colini

	vmin = fiso(1)
	vmax = fiso(isoanz)

	end

c************************************************************

	subroutine colentry(icol,viso,col)

c returns entry of color table

	use color

	implicit none

	integer icol
	real viso,col

	call colini

	viso = 0.
	col = 0.

	if( icol .le. 0 ) then
	  return
	else if( icol .le. isoanz ) then
	  viso = fiso(icol)
	  col = ciso(icol)
	else if( icol .eq. isoanz+1 ) then
	  col = ciso(icol)
	else
	  return
	end if

	end
 
c************************************************************

	function getcol(value)

c gets color for value

	use color

	implicit none

	real getcol
	real value
	integer i

	do i=1,isoanz
	  getcol = ciso(i)
	  if( fiso(i) .gt. value ) return
	end do

	getcol = ciso(isoanz+1)

	end

c************************************************************

	function getcolval(rindex)

c gets value for real index
c
c if array has been read -> closest value
c if regular values (valmin/max) -> use real rindex

	use color

	implicit none

	real getcolval
	real rindex

	real ri,val,dval

	ri = rindex
	ri = min(ri,float(isoanz))
	ri = max(ri,1.)

	if( nisord .gt. 0 .or. isoanz .eq. 1 ) then	!values read
	  val = fiso(nint(ri))
	else if( iusear .ne. 0 ) then			!use array
	  val = fiso(nint(ri))
	else
	  dval = (fiso(isoanz) - fiso(1)) / (isoanz-1)
	  val = fiso(1) + (ri-1.) * dval
	end if

	getcolval = val

	end
 
c************************************************************

	subroutine make_single_isolines(ntick)

c computes values for which single isolines are plotted

	use color

	implicit none

	integer ntick

	integer i
	real dtick,rit,value
	real getcolval

	if( ntick .le. 1 ) then
	  nriso = isoanz
	  do i=1,nriso
	    riso(i) = fiso(i)
	  end do
	else
	  nriso = ntick
	  dtick = float(isoanz-1) / (nriso-1)
	  do i=1,nriso
            rit = 1. + (i-1)*dtick
            value = getcolval(rit)
	    riso(i) = value
	  end do
	end if

	write(6,*) 'make_single_isolines: ',nriso
	write(6,*) (riso(i),i=1,nriso)

	end

c************************************************************
c************************************************************
c************************************************************

	subroutine set_default_color_table( ictab )

c sets default color table

	implicit none

	integer ictab

	call qtdef(ictab)

	end

c*****************************************************************

	subroutine get_color_table( icolor )

c gets color table used

	use color

	implicit none

	integer icolor

	icolor = icol_def

	end

c*****************************************************************

	subroutine set_color_table( icolor )

c sets color table to be used

	use color

	implicit none

	integer icolor

	icol_def = icolor
	call qsetc(0.) 		!try if color table exists

	end

c*****************************************************************

	subroutine reset_color_table

c sets default color table

	use color

	implicit none

	call set_color_table( icol_def )

	end

c*****************************************************************

	subroutine write_color_table

c writes default color table to terminal

	use color

	implicit none

	write(6,*) 'color table: actual = ',icol_def
	!write(666,*) 'color table: actual = ',icol_def

	end

c*****************************************************************
c
c set image colorscale hsb 0.666 0.0 1.0 .min. hsb 0.666 1.0 1.0 .max.
c Panel "3. White-blue (HSB blending)"
c 
c set image colorscale rgb 1.0 0.0 0.0 .min. rgb 0.0 0.0 1.0 .max.
c Panel "4a. Red-blue (RGB blending)"
c
c*****************************************************************

	subroutine qsetc( col )

c changes color using the actual color table

	use color

	implicit none

	real col

	integer icol

	icol = icol_def

	if( icol .eq. -1 ) then
	  call qcolor(col)
	else if( icol .eq. 0 ) then
	  call qgray(col)
	else if( icol .eq. 1 ) then
	  call qhue(col)
	else if( icol .eq. 2 ) then
	  call white_blue( col )
	else if( icol .eq. 3 ) then
	  call white_red( col )
	else if( icol .eq. 4 ) then
	  call blue_white_red( col )
	else if( icol .eq. 5 ) then
	  call blue_black_red( col )
	else if( icol .eq. 6 ) then
	  call hue_sqrt( col )
	else if( icol .eq. 7 ) then
	  call hue_pow2( col )
	else if( icol .eq. 8 ) then
	  call col_cust( col )
	else
	  write(6,*) 'icol = ',icol
	  stop 'error stop qsetc: no such color table'
	end if

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine white_blue( color )
	implicit none
	real color
	call qhsb(0.666,color,1.)
	end

	subroutine white_red( color )
	implicit none
	real color
	call qhsb(1.,color,1.)
	end

	subroutine blue_white_red( color )
	implicit none
	real color
	if( color .le. 0.5 ) then
	  call qhsb(0.666,1.-2.*color,1.)
	else
	  call qhsb(1.,2.*(color-0.5),1.)
	end if
	end

	subroutine blue_black_red( color )
	implicit none
	real color
	if( color .le. 0.5 ) then
	  call qhsb(0.666,1.,1.-2.*color)
	else
	  call qhsb(1.,1.,2.*(color-0.5))
	end if
	end

	subroutine hue_sqrt( color )
	implicit none
	real color
	call qhue(sqrt(color))
	end

	subroutine hue_pow2( color )
	implicit none
	real color
	call qhue(color*color)
	end

	subroutine col_cust( col )
	use color
	implicit none
	real col
	real c
	integer ic
	c = col
	c = min(1.,c)
	c = max(0.,c)
	ic = nint((icmax-1)*c) + 1
	call qrgb(coltab(1,ic),coltab(2,ic),coltab(3,ic))
	!write(6,*) 'color... : ',color,c,ic,coltab(:,ic)
	end

c*****************************************************************

