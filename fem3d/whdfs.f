c
c $Id: whdfs.f,v 1.2 1998/01/22 16:50:15 georg Exp $

	program whdf

c writes hdf files

	parameter (irank=3)
	parameter (idimx=25,idimy=9,idimz=200)

	real s(idimx,idimy,idimz)

	character*80 line
	integer iaux(25)

	lu=0
	lv=0
	lw=0
	ls=0

	write(6,*) 'This program reads from stdin'

    1	continue

	read(5,'(a)',end=99) line

	if(line(1:4).eq.'salt') then

	  ls=ls+1
	  if(ls.gt.idimz) stop 'error dim z'
	  do iy=1,idimy
	    read(5,'(a)') line
	    read(line,'(i2,2x,25i3)') idum,iaux
	    do ix=1,idimx
	      s(ix,iy,ls) = iaux(ix)
	    end do
	  end do

        else 

	  do iy=1,idimy
	    read(5,'(a)') line
	  end do

	end if

	goto 1

   99	continue

	write(6,*) lu,lv,lw,ls

	iaux(1)=idimx
	iaux(2)=idimy
	iaux(3)=ls
	iret = dsadata('sss.hdf',irank,iaux,s)

	stop
	end
