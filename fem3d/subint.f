
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

	subroutine intsha(x,y,z,icolor,rlev,ncol,itrat)
c
c interpolation in triangle and area plot of colors
c
	dimension x(3),y(3),z(3)
	dimension icolor(1),rlev(1)
c
	logical btypa
	dimension xp(7),yp(7)
	dimension xph(7),yph(7) !only because qareac changes vector
				!-->change  qareac
c
	save eps
	data eps /1.e-5/
c
c find min/max %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	call mima(z,3,zmin,zmax)
c
c put nodes in order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c a is first node, b is last node
c 1 is first free node, 2 (only if btypa=.false.) is second free node
c
c find first and last node
c
	icont=0
	ieck=0
	do ii=3,1,-1
	if(z(ii).eq.zmin) then
	    if(icont.eq.0) then
		xa=x(ii)
		ya=y(ii)
		za=z(ii)
		icont=icont+1
		ieck=ieck+ii
	    else if(icont.eq.1) then
		xb=x(ii)
		yb=y(ii)
		zb=z(ii)
		icont=icont+1
		ieck=ieck+ii
	    else if(icont.eq.2) then
		icont=icont+1
	    end if
	end if
	end do
c
c find free nodes
c
	if(icont.eq.1) then
		xb=xa
		yb=ya
		zb=za
		ii=mod(ieck,3)+1
		x1=x(ii)
		y1=y(ii)
		z1=z(ii)
		ii=mod(ii,3)+1
		x2=x(ii)
		y2=y(ii)
		z2=z(ii)
		btypa=.false.
	else if(icont.eq.2) then
		ii=6-ieck
		x1=x(ii)
		y1=y(ii)
		z1=z(ii)
		btypa=.true.
	else if(icont.eq.3) then
		ii=1
		x1=x(ii)
		y1=y(ii)
		z1=z(ii)
		btypa=.true.
	end if
c
c find first level .ge. nodes a and b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do icol=1,ncol-1
	if(rlev(icol).ge.zmin) goto 2
	end do
    2   continue
c
c start plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	do while(icol.le.ncol)
c
c zr is actual level (last level must be always true)
c
	if(icol.lt.ncol) then
		zr=rlev(icol)
	else
		zr=zmax+1.
	end if
c
c two cases
c
	if(btypa) then
c
c case a (one free node)
c
		if(z1.le.zr) then
			xp(1)=xa
			yp(1)=ya
			xp(2)=x1
			yp(2)=y1
			xp(3)=xb
			yp(3)=yb
			ip=3
			call qareac(xp,yp,ip,icolor(icol),itrat)
			goto 99
		else
			xp(1)=xa
			yp(1)=ya
			dz=z1-za
			if(dz.gt.eps) then
				xp(2)=xa+(x1-xa)*(zr-za)/dz
				yp(2)=ya+(y1-ya)*(zr-za)/dz
			else
				xp(2)=(xa+x1)*0.5
				yp(2)=(ya+y1)*0.5
			end if
			dz=z1-zb
			if(dz.gt.eps) then
				xp(3)=xb+(x1-xb)*(zr-zb)/dz
				yp(3)=yb+(y1-yb)*(zr-zb)/dz
			else
				xp(3)=(xb+x1)*0.5
				yp(3)=(yb+y1)*0.5
			end if
			xp(4)=xb
			yp(4)=yb
			ip=4
			do i=1,ip
			xph(i)=xp(i)
			yph(i)=yp(i)
			end do
			call qareac(xp,yp,ip,icolor(icol),itrat)
			xa=xph(2)
			ya=yph(2)
			za=zr
			xb=xph(3)
			yb=yph(3)
			zb=zr
		end if
	else
c
c case b (two free nodes)
c
		xp(1)=xa
		yp(1)=ya
		ip=1
		ieck=0
		if(z1.gt.zr) then
			ip=ip+1
			dz=z1-za
			if(dz.gt.eps) then
				xp(ip)=xa+(x1-xa)*(zr-za)/dz
				yp(ip)=ya+(y1-ya)*(zr-za)/dz
			else
				xp(ip)=(xa+x1)*0.5
				yp(ip)=(ya+y1)*0.5
			end if
			ieck=ieck+1
			ipa=ip  !ipa is new a
		else
			ip=ip+1
			xp(ip)=x1
			yp(ip)=y1
		end if
		if(z2.gt.zr.and.z1.le.zr) then
			ip=ip+1
			dz=z2-z1
			if(dz.gt.eps) then
				xp(ip)=x1+(x2-x1)*(zr-z1)/dz
				yp(ip)=y1+(y2-y1)*(zr-z1)/dz
			else
				xp(ip)=(x2+x1)*0.5
				yp(ip)=(y2+y1)*0.5
			end if
			ipf=2   !ipf is new free node
			ipa=ip  !ipa is new a
		else if(z1.gt.zr.and.z2.le.zr) then
			ip=ip+1
			dz=z2-z1
			if(-dz.gt.eps) then
				xp(ip)=x1+(x2-x1)*(zr-z1)/dz
				yp(ip)=y1+(y2-y1)*(zr-z1)/dz
			else
				xp(ip)=(x2+x1)*0.5
				yp(ip)=(y2+y1)*0.5
			end if
			ipf=1   !ipf is new free node
			ipb=ip  !ipb is new b
		else if(z2.le.zr) then
			ip=ip+1
			xp(ip)=x2
			yp(ip)=y2
		end if
		if(z2.gt.zr) then
			ip=ip+1
			dz=z2-zb
			if(dz.gt.eps) then
				xp(ip)=xb+(x2-xb)*(zr-zb)/dz
				yp(ip)=yb+(y2-yb)*(zr-zb)/dz
			else
				xp(ip)=(xb+x2)*0.5
				yp(ip)=(yb+y2)*0.5
			end if
			ieck=ieck+1
			ipb=ip  !ipb is new b
		else
			ip=ip+1
			xp(ip)=x2
			yp(ip)=y2
		end if
		ip=ip+1
		xp(ip)=xb
		yp(ip)=yb
c
		do i=1,ip
		xph(i)=xp(i)
		yph(i)=yp(i)
		end do
c
		call qareac(xp,yp,ip,icolor(icol),itrat)
c
		if(ieck.eq.0) then
			goto 99
		else
			xa=xph(ipa)
			ya=yph(ipa)
			za=zr
			xb=xph(ipb)
			yb=yph(ipb)
			zb=zr
			if(ieck.ne.2) then
				btypa=.true.
				if(ipf.eq.2) then
					x1=x2
					y1=y2
					z1=z2
				end if
			end if
		end if
	end if
c
	icol=icol+1
c
	end do
c
   99   continue
c
	return
	end
