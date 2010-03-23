
c subroutine axis(xpage,ypage,ibcd,nchar,axlen,angle,firstv,deltav)
c subroutine scale(array,axlen,npts,inc)
c subroutine line(xarray,yarray,npts,inc,lintyp,inteq)
c subroutine number(xpage,ypage,height,fpn,angle,ndec)

c ggu 03.04.2009 CW -> avoid compiler warning

      subroutine axis(xpage,ypage,ibcd,nchar,axlen,angle,firstv,deltav)

	dimension ibcd(1)
	!data isym/4h*10 /

      character symb*40
      integer isym(40)          !ggu 03.04.2009 CW -> avoid compiler warning
      equivalence (symb,isym(1))

      do k=1,10
        symb(k:k)=' '
      end do
c
      kn=nchar
      a=1.0
c set constants for annotation on cw or ccw side of axis
      if (kn) 1,2,2
1     a=-a
      kn=-kn
2     ex=0.0
      adx= abs  (deltav)
      if (adx) 3,7,3
3     if (adx- 99.0) 6,4,4
4     adx=adx/10.0
      ex=ex+1.0
      go to 3
5     adx=adx*10.0
      ex=ex-1.0
6     if (adx-0.01) 5,7,7
7     xval=firstv*10.0**(-ex)
      adx=2.0*(deltav*10.0**(-ex))
      sth=angle*0.0174533
      cth=cos(sth)
      sth=sin(sth)
      cth2=cth+cth
      sth2=sth+sth
      dxb=-0.254
      dyb=0.38*a-0.127
      xn=xpage+dxb*cth-dyb*sth
      yn=ypage+dyb*cth+dxb*sth
      ntic=axlen+1.0
      ntc=(ntic+1)/2
      nt=ntc/2
      do 20  i=1,ntc
      call number(xn,yn,0.24,xval,angle,2)
      if(i-1) 10,10,105
10    xn=xpage+2.0*dxb*cth-dyb*sth
      yn=ypage+2.0*dxb*sth+dyb*cth
105   xval=xval+adx
      xn=xn+cth2
      yn=yn+sth2
      if (nt) 20,11,20
11    z=kn
      if (ex)  12,13,12
12    z=z+7.0
13    dxb=-.175*z+axlen*0.5
      dyb=0.8*a-0.2
      xt=xpage+dxb*cth-dyb*sth
      yt=ypage+dyb*cth+dxb*sth
      call symbol(xt,yt,0.36,ibcd(1),angle,kn)
      if (ex)  14,20,14
14    z=kn+2
      xt=xt+z*cth*0.36
      yt=yt+z*sth*0.36
      call symbol(xt,yt,0.36,isym,angle,3)
      xt=xt+(3.0*cth-0.6*sth)*0.36
      yt=yt+(3.0*sth+0.6*cth)*0.36
      call number(xt,yt,0.18,ex,angle,-1)
20    nt=nt-1
      call plot(xpage+axlen*cth,ypage+axlen*sth,3)
      dxb=-0.178*a*sth
      dyb=+0.178*a*cth
      a=ntic-1
      xn=xpage+a*cth
      yn=ypage+a*sth
      do  30  i=1 , ntic
      call plot(xn,yn,2)
      call plot(xn+dxb,yn+dyb,2)
      call plot(xn,yn,2)
      xn=xn-cth
      yn=yn-sth
30    continue
      return
      end




      subroutine scale(array,axlen,npts,inc)

      dimension  array(1),save(7)

      save(1)= 1.0
      save(2)= 2.0
      save(3)= 4.0
      save(4)= 5.0
      save(5)= 8.0
      save(6)=10.0
      save(7)=20.0
      fad=0.01
      k=iabs(inc)
      n=npts*k
      y0=array(1)
      yn=y0
      do  25  i=1 ,n,k
      ys=array(i)
      if  (y0-ys)  22,22,21
21    y0=ys
      go  to  25
22    if  (ys-yn)  25,25,24
24    yn=ys
25    continue
      firstv=y0
      if  (y0)  34,35,35
34    fad=fad-1.0
35    deltav=(yn-firstv)/axlen
      if (deltav) 70,70,40
40    i=alog10(deltav) + 1000.0
      p=10.0**(i-1000)
      deltav=deltav/p-0.01
      do  45  i=1,6
      is=i
      if  (save(i)-deltav)  45,50,50
45    continue
50    deltav=save(is)*p
      firstv=deltav*aint(y0/deltav+fad)
      t=firstv+(axlen+0.01)*deltav
      if (t-yn)  55,57,57
55    firstv=p*aint(y0/p+fad)
      t=firstv+(axlen+.01)*deltav
      if (t-yn)  56,57,57
56    is=is+1
      go  to  50
57    firstv=firstv-aint((axlen+(firstv-yn)/deltav)/2.0)*deltav
      if (y0*firstv) 58,58,59
58    firstv=0.0
59    if  (inc) 61,61,65
61    firstv=firstv+aint(axlen+.5)*deltav
      deltav=-deltav
65    n=n+1
      array(n)=firstv
      n=n+k
      array(n)=deltav
67    return
70    deltav=2.0*firstv
      deltav=abs(deltav/axlen)+1.
      go to 40

      end




      subroutine line(xarray,yarray,npts,inc,lintyp,inteq)

      dimension xarray(1),yarray(1)

      lmin = npts*inc+1
      ldx  = lmin+inc
      nl   = lmin-inc
      firstx = xarray(lmin)
      deltax = xarray(ldx)
      firsty = yarray(lmin)
      deltay = yarray(ldx)
      call where (xn,yn,df)
      df=amax1(abs((xarray( 1)-firstx)/deltax-xn),
     1         abs((yarray( 1)-firsty)/deltay-yn) )
      dl=amax1(abs((xarray(nl)-firstx)/deltax-xn),
     1         abs((yarray(nl)-firsty)/deltay-yn) )

      ipen = 3
      icode = -1
      nt =iabs(lintyp)

      if(lintyp.eq.0) nt = 1
      if(df.gt.dl) then
        nf = nl
        na = ((npts-1)/nt)*nt+nt-(npts-1)
        kk = -inc
      else
        nf = 1
        na = nt
        kk = inc
      end if

      if(lintyp.lt.0) then
        ipena = 3
        icodea = -1
        lsw = 1
      else
        if(lintyp.eq.0) na = ldx
        ipena = 2
        icodea = -2
        lsw=0
      end if

      do i = 1,npts
        xn = (xarray(nf)-firstx)/deltax
        yn = (yarray(nf)-firsty)/deltay
        if(na.eq.nt) then
          call symbol (xn,yn,0.20,inteq,0.0,icode)
          na = 1
        else if(na.lt.nt) then
          if(lsw.eq.0) call plot (xn,yn,ipen)
          na = na + 1
        else
          call plot (xn,yn,ipen)
          na = na + 1
        end if
        nf = nf+kk
        icode = icodea
        ipen = ipena
      end do

      return
      end



      subroutine number(xpage,ypage,height,fpn,angle,ndec)

      character num*20,iar*10,minus*1,ipoint*1,blank*1
      integer inum(20)          !ggu 30.11.2001 CW -> avoid compiler warning
      equivalence (num,inum(1)) !ggu 30.11.2001 CW -> avoid compiler warning
      data iar/'0123456789'/
      data minus/'-'/,ipoint/'.'/,blank/' '/

      do k=1,20
        num(k:k)=blank
      end do

      ii=0
      fpv = fpn
      n = ndec
      maxn = 9

      if(n.gt.maxn) n=maxn
      if(-n.gt.maxn) n=-maxn

      if(fpv.lt.0) then
        ii=ii+1
        num(ii:ii)=minus
      end if

      mn = -n
      if(n.lt.0) mn = mn - 1
      fpv = abs(fpv) + (0.5 * 10. ** mn)
      i = alog10(fpv) + 1.0

      ilp = i
      if(n.lt.-1) ilp = ilp + n + 1

      if(ilp.le.0) then
        ii=ii+1
        num(ii:ii)=iar(1:1)
      else
        if(ilp+n.gt.18) then
          n = -1
          if(ilp.gt.19) ilp = 19
        end if

        do j=1,ilp
          k = fpv * 10. ** (j - i)
          if(k.gt.9) k = 9
          ii=ii+1
          isb=k+1
          num(ii:ii)=iar(isb:isb)
          fpv = fpv - (float(k) * 10. ** (i - j))
        end do
      end if

      if(n.ge.0) then
        ii=ii+1
        num(ii:ii)=ipoint
        if(n.gt.0) then
          do j=1,n
            k = fpv * 10.
            if(k.gt.9) k = 9
            ii=ii+1
            isb=k+1
            num(ii:ii)=iar(isb:isb)
            fpv = fpv * 10. - float(k)
          end do
        end if
      end if

c     nchar=ii+1000
      nchar=ii

      !call symbol (xpage,ypage,height,num,angle,nchar)
      call symbol (xpage,ypage,height,inum,angle,nchar) !!ggu 30.11.2001 CW

      return
      end
