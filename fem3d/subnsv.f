c
c $Id: subnsv.f,v 1.3 2000/05/26 11:22:08 georg Exp $
c
c routines for pre processing routines
c
c contents :
c
c subroutine sp190(nen,nel,nform,mbw)           determine bandwidth mbw
c
c revision log :
c
c 15.05.2000    ggu     subroutine sp191 removed
c
c**********************************************************

	subroutine sp190(nen,nel,nform,mbw)

c determine bandwidth mbw
c
c nen           node index
c nel           number of elements
c nform         number of nodes per element
c mbw           bandwidth (return value)

	dimension nen(1)

	mh=0

	do ie=1,nel
	do i=1,nform
	k=nen(nform*(ie-1)+i)
	if(k.gt.0) then
		do ii=i+1,nform
		kk=nen(nform*(ie-1)+ii)
		if(kk.gt.0) then
			mm=iabs(kk-k)
			if(mm.gt.mh) mh=mm
		end if
		end do
	end if
	end do
	end do

	mbw=mh

	end

c**********************************************************

