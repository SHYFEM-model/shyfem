c
c $Id: whdf1.f,v 1.2 1998/01/22 16:50:14 georg Exp $

	program whdf

c writes hdf files

	parameter (irank=3)

	real a(25,6,20)
	integer idim(irank)
	data idim /25,6,0/

    1	continue
	read(5,*) it
	if( it.lt.9000 ) goto 1

	l=0
    2	continue
	l=l+1
	do i=6,1,-1
	  read(5,*) iaux,(a(j,i,l),j=1,25)
	end do
	read(5,*,end=99) it
	if(l.lt.10) goto 2

   99	continue

	idim(3)=l
	iret = dsadata('example.hdf',irank,idim,a)

	stop
	end
