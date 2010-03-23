c
c $Id: subfnm.f,v 1.7 1998/06/23 09:09:21 georg Exp $
c
c fnm text retrieving routines
c
c contents :
c
c function manfnm(name,text,id,what)    internal routine to access fnm data
c subroutine getfnm(name,text)          gets string name and writes it in text
c function itsfnm(name)                 tests if name is available
c function iscfnm(name,sect)		tests if name is in section 
c subroutine sctfnm(sect)		sets default name of section
c subroutine putfnm(name,text)          puts string text in name
c subroutine addfnm(name,text)          adds string text to name
c function inffnm(type)                 info about parameters
c subroutine lstfnm(name,text,id)       lists entry id in name and text
c subroutine chafnm                     changes string interactively
c subroutine prifnm                     prints parameter values
c
c revision log :
c
c revised 04.04.92 by ggu (cint)
c revised 12.06.97 by ggu (section introduced)
c 18.03.1998	ggu	save secfnm (bug uncovered by g77)
c
c**************************************************************
c
	function manfnm(name,text,id,what)
c
c internal routine to access fnm data
c
c name          name of string
c text          contents of string
c id            number of entry (only for what='l')
c what          what to do ?
c               g : return contents of name in text
c               t : name is defined ?
c               s : put section of name in auxsec
c               p : put text into name
c               a : add text with name
c               f : how many entries ?
c               l : list name and text of entry id
c manfnm        entry number (0 if no entry name is found,
c               ...total number of entries for what='f')
c
c comments:
c
c - uses value in parsec for section of name
c
c problems:
c
c - there is no way to change section for an already given parameter
c
c
	character*(*) name,text,what
c-------------------------------------------
	parameter (ntxtdi=1000,nnamdi=100)
c-------------------------------------------
	character*1 txtfnm(ntxtdi)      !array containing text strings
	integer iptfnm(nnamdi+1)        !array with pointers to strings
	character*6 namfnm(nnamdi)      !array containing names of strings
        character*6 secfnm(nnamdi)      !array containing section names
c-------------------------------------------
        character*6 fnmsec, auxsec
        common /fnmsco/ fnmsec, auxsec          !default section
c-------------------------------------------
c
	character*6 namarg,namvec,actsec
	character*1 type,cint
c	logical bdebug
c
	save ifirst,nentry
	save txtfnm,iptfnm,namfnm,secfnm
	data ifirst,nentry / 0 , 0 /
c
c	bdebug = .false.
c
	if(ifirst.eq.0) then    !initialize pointer
		fnmsec = ' '
		iptfnm(1)=1
		ifirst=1
	end if
c
	type=what(1:1)
	call uplow(type,'low')
c
	namarg=name
	call uplow(namarg,'low')
c
c find entry name %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	ientry=0
	do i=1,nentry
	   namvec=namfnm(i)
	   call uplow(namvec,'low')
	   if(namarg.eq.namvec) ientry=i
	end do
c
c see what to do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(type.eq.'f') then                    !filling of array
		manfnm=nentry
		return
	else if(type.eq.'t') then               !test if name exists
		manfnm=ientry
		return
	else if(type.eq.'s') then               !put section of name in auxsec
		auxsec = ' '
		if( ientry .gt. 0 ) then
		  auxsec=secfnm(ientry)
		end if
		manfnm=ientry
	        return
	else if(type.eq.'l') then               !list id entry
		text=' '
		if(id.le.0.or.id.gt.nentry) then
			manfnm=0
			name=' '
			return
		else
			ientry=id
			manfnm=ientry
			name=namfnm(ientry)
		end if
	else if(type.eq.'g') then               !get string name
		text=' '
		if(ientry.eq.0) then
			manfnm=0
			return
		else
			manfnm=ientry
		end if
	else if(type.eq.'a') then               !add name
		if(ientry.eq.0) then
			manfnm=nentry+1
		else
			manfnm=ientry
		end if
	else if(type.eq.'p') then               !put name
		if(ientry.eq.0) then
			manfnm=0
			return
		else
			manfnm=ientry
		end if
	else                                    !not recognized
		write(6,*) 'type not recognized : ',what
		stop 'error stop : manfnm'
	end if
c
c read or write text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	if(type.eq.'g'.or.type.eq.'l') then             !get text
		ianf=iptfnm(ientry)
		iend=iptfnm(ientry+1)-1
		ii=0
		do i=ianf,iend
		   ii=ii+1
		   cint=txtfnm(i)
		   text(ii:ii)=cint
c		   text(ii:ii)=char(int)
		end do
	else if(type.eq.'a'.or.type.eq.'p') then        !put text
		ianz=ichanm(text)
		actsec = fnmsec
                if(ientry.ne.0) then                    !see if too big
			ianf=iptfnm(ientry)
			iend=iptfnm(ientry+1)-1
			actsec = secfnm(ientry)		!save for later
			if(ianz.gt.iend-ianf+1) then    !too big
				namfnm(ientry)=' '
				secfnm(ientry)=' '
				do i=ianf,iend
				   txtfnm(i)=' '
				end do
				ientry=0
				manfnm=nentry+1
			end if
		end if
		if(ientry.eq.0) then                    !new string -> at end
			nentry=nentry+1
			ientry=nentry
			if(ientry.gt.nnamdi) then
				write(6,*) 'dimension error nnamdi'
				stop 'error stop : manfnm'
			end if
			ianf=iptfnm(ientry)
			iend=ianf+ianz-1
			if(iend.gt.ntxtdi) then
				write(6,*) 'dimension error ntxtdi'
				stop 'error stop : manfnm'
			end if
			namfnm(ientry)=name
			secfnm(ientry) = actsec
			iptfnm(ientry+1)=iend+1
			ii=0
			do i=ianf,iend
			ii=ii+1
			cint=text(ii:ii)
			txtfnm(i)=cint
c			int=ichar(text(ii:ii))
c			txtfnm(i)=int
			end do
		else                                    !put in old place
			iact=ianf
			do ii=1,ianz
			cint=text(ii:ii)
			txtfnm(iact)=cint
c			int=ichar(text(ii:ii))
c			txtfnm(iact)=int
			iact=iact+1
			end do
			do i=iact,iend
			txtfnm(i)=' '
			end do
		end if
	end if
c
	return
	end
c
c***********************************************************
c
	subroutine getfnm(name,text)
c
c gets string name and writes it in text
c
	character*(*) name,text
c
	i=manfnm(name,text,0,'get')
	if(i.eq.0) then
		write(6,*) 'String not found : ',name
		stop 'error stop : getfnm'
	end if
c
	return
	end
c
c***********************************************************
c
	function itsfnm(name)
c
c tests if name is available (itsfnm =|= 0) or not (itsfnm=0)
c
	character*(*) name
	character*6 text
c
	itsfnm=manfnm(name,text,0,'test')
c
	return
	end
c
c***********************************************************
c
	function iscfnm(name,sect)
c
c tests if name is in section (iscfnm =|= 0) or not (iscfnm=0)
c
	character*(*) name,sect
	character*6 text,asect
c-------------------------------------------
	character*6 fnmsec, auxsec
	common /fnmsco/ fnmsec, auxsec		!default section
c-------------------------------------------
c
c	next call puts section of name into auxsec
c
	iscfnm=manfnm(name,text,0,'sect')
c
	if( iscfnm .gt. 0 ) then
	    asect = sect
	    call uplow(asect,'low')
	    if( asect .ne. auxsec ) iscfnm = 0
	end if
c
	return
	end
c
c***********************************************************
c
	subroutine sctfnm(sect)
c
c sets default name of section
c
	character*(*) sect
	character*6 name,text
c-------------------------------------------
	character*6 fnmsec, auxsec
	common /fnmsco/ fnmsec, auxsec		!default section
c-------------------------------------------
c
c	next call is dummy call to be sure that fnmsco is initialized
c
	idum=manfnm(name,text,0,'fill')
	idum = 2 * idum
c
	auxsec = sect
	call uplow(auxsec,'low')
	fnmsec = auxsec
c	write(6,*) 'sctfnm: default section = ',fnmsec
c
	return
	end
c
c***********************************************************
c
	subroutine putfnm(name,text)
c
c puts string text in name
c
	character*(*) name,text
c
	i=manfnm(name,text,0,'put')
	if(i.eq.0) then
		write(6,*) 'String not found : ',name
		stop 'error stop : putfnm'
	end if
c
	return
	end
c
c**************************************************************
c
	subroutine addfnm(name,text)
c
c adds string text to name
c
	character*(*) name,text
c
	idum=manfnm(name,text,0,'add')
	idum=idum*2
c
	return
	end
c
c**************************************************************
c
	function inffnm(type)
c
c info about parameters :
c       type = 'f' --> how many entries
c
	character*(*) type
	character*6 name,text
c
	if(type(1:1).eq.'f') then       !how many entries
		inffnm=manfnm(name,text,0,'fill')
	else
		write(6,*) 'type not recognized : ',type
		stop 'error stop : inffnm'
	end if
c
	return
	end
c
c**************************************************************
c
	subroutine lstfnm(name,text,id)
c
c lists entry id in name and text
c
	character*(*) name,text
c
	idum=manfnm(name,text,id,'list')
	idum=idum*2
c
	return
	end
c
c*****************************************************
c
	subroutine chafnm
c
c changes string interactively
c
	implicit none
c
c local
	logical lparam
	character*6 name
	character*80 text,text1
	integer nlen
c function
	integer itsfnm,ichanm
c
	lparam=.true.
c
	do while(lparam)
		call prifnm
c
		write(6,*) 'Enter name of string :'
		read(5,'(a)') name
c
		call uplow(name,'low')
c
		if(name.ne.' ') then
			if(itsfnm(name).eq.0) then
				write(6,*)'No string with this name'
			else
			   call getfnm(name,text)
			   nlen=ichanm(text)
			   if(nlen.le.0) nlen=1
			   write(6,'(1x,a6,a1,a)') name,'=',text(1:nlen)
			   write(6,*) 'Enter new string : '
     +                          ,'(<CR> for OK, - to delete)'
			   read(5,'(a)') text1
			   if(ichanm(text1).le.0) then
				text1=text
			   else if(text1(1:1).eq.'-') then
				text1=' '
			   end if
			   call putfnm(name,text1)
			end if
		else
			lparam=.false.
		end if
c
	end do
c
	return
	end

c*****************************************************

	subroutine prifnm

c prints parameter values

	implicit none

c common
c-------------------------------------------
	character*6 fnmsec, auxsec
	common /fnmsco/ fnmsec, auxsec		!default section
c-------------------------------------------
c local
	character*80 text
	character*6 name
	integer i,nentry,nlen
c function
	integer inffnm,ichanm
	integer manfnm

	nentry=inffnm('fill')

	do i=1,nentry
	  call lstfnm(name,text,i)
          if(name.ne.' ') then
	    nlen=manfnm(name,text,0,'sect')	!nlen is dummy
	    nlen=ichanm(text)
	    if(nlen.le.0) nlen=1
	    write(6,2345) i,nlen,name,auxsec,text(1:nlen)
          end if
	end do

	return
 2345   format(1x,2i4,2(1x,a6,1x),3x,a)
	end
