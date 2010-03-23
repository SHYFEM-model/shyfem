c
c $Id: wincrea.f,v 1.4 1998/06/18 10:39:25 georg Exp $
c
        program wincrea

c This program contains an example of how to create a wind file
c to be used with SHYFEM. The wind is always constant in space.
c If you need to create a spatially non constant wind please
c refer to the documentation.
c
c To create a wind file constant in time only one wind record
c has to be created. In this case the time stamp for the wind
c record is meaningless. In the case of wind variable in time
c more records have to be created that describe the time
c evolution of the wind. The time stamps of the wind records
c refer to the time used by the simulation.
c
c There is no need to create a wind record for every time
c step in the simulation. It is enough that the wind records
c describe faithfully the evolution of the wind. Between the
c given time records the wind is interpolated linearly by
c the simulation. 
c
c However, be sure that there are enough time records
c to cover the whole period of the simulation.
c
c To be used with SHYFEM, please insert the following lines in the
c in the STR file :
c
c $name
c 	...
c 	wind='/usr/users/model/wind/scir.win'
c 	...
c $end
c
c in section $name for the name of the wind file (change the 
c directory and file name to suit them to your case) and
c
c $para
c 	...
c	dragco = 3.2e-3
c 	...
c $end
c
c in section $para for the value of the drag coefficient.
c Drag coefficients between 1.8e-3 and 3.2e-3 are realistic.
c
c A call to "mkwind" writes one wind record to the wind file.
c The call to mkwind is as follows
c
c	call mkwind(iunit,itime,speed,dir)
c
c with the following meaning
c
c iunit		unit for output
c itime		record refers to time itime (in seconds)
c speed		wind speed in m/s [0 - ...]
c dir		wind direction [0 - 360]
c		... 0. from north, 90. from east...
c
c examples for dir :
c
c bora:		dir = 45.
c scirocco:	dir = 135.
c
c The following example shows how to create a wind file
c for scirocco wind of 10 m/s constant in time.
c
c	open(1,file='scir.win',form='unformatted',status='new')
c	call mkwind(1,0,10.,135.)
c	close(1)
c
c The next example shows how to create a wind file
c for bora wind. The wind speed is 15 m/s at time 0
c and falls off to 5 m/s at time 43200 (12 hours)
c to stay constant for the rest of the day.
c
c	open(1,file='bora.win',form='unformatted',status='new')
c	call mkwind(1,0,15.,45.)
c	call mkwind(1,43200,5.,45.)
c	call mkwind(1,86400,5.,45.)
c	close(1)
c
c The next example shows a bora wind (15 m/s) that
c gradually during 12 hours turns into a scirocco wind 
c and falls off to 7 m/s.
c
c	open(1,file='mixed.win',form='unformatted',status='new')
c	call mkwind(1,0,15.,45.)
c	call mkwind(1,21600,11.,90.)
c	call mkwind(1,43200,7.,135.)
c	close(1)
c
c	open(1,file='mixed.win',form='unformatted',status='new')
c	call mkwind(1,0,15.,45.)
c	call mkwind(1,21600,11.,90.)
c	call mkwind(1,43200,7.,135.)
c	close(1)
c
c	open(1,file='south.win',form='unformatted',status='unknown')
c	call mkwind(1,0,10.,180.)
c	close(1)

	open(1,file='bora.win',form='unformatted',status='unknown')
	call mkwind(1,0,10.,45.)
	close(1)

	open(1,file='test.win',form='unformatted',status='unknown')
	call mkwind(1,0,20.,0.)
	close(1)

        end

c********************************************************************

	subroutine mkwind(iunit,itime,speed,dir)

c writes one record of wind data
c
c iunit		unit for output
c itime		record refers to time itime (in seconds)
c speed		wind speed in m/s [0 - ...]
c dir		wind direction [0 - 360]
c		... 0. from north, 90. from east...

	implicit none

	integer iunit,itime
	real speed,dir

	integer one
	parameter(one=1)

	real wx,wy

	call cvwind(speed,dir,wx,wy)

        write(iunit) itime,one
        write(iunit) wx,wy

	return
	end

c********************************************************************

	subroutine cvwind(speed,dir,wx,wy)

c converts polar wind data into cartesian format
c
c speed		wind speed in m/s [0 - ...]
c dir		wind direction [0 - 360]
c		... 0. from north, 90. from east...
c wx,wy		wind vector in cartesian coordinates (return)

	implicit none

	real speed,dir
	real wx,wy

c alphla is angle between true north and y-coordinate of lagoon        
c ... 0. for venlag61, 90. for bas059, 74.3 for lip

	real alphla
	real pi,rad
	parameter(alphla=0.)
	parameter(pi=3.14159,rad=pi/180.)

	real gamma

        gamma = rad*( 270. - alphla - dir )
        wx=speed*cos(gamma)
        wy=speed*sin(gamma)

	return
	end
