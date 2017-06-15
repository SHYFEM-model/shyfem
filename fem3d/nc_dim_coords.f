
c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine get_nc_dimensions(ncid,bverb,nt,nx,ny,nz)

        implicit none

        integer ncid
        logical bverb
        integer nt,nx,ny,nz

        integer dim_id,n,ndims,time_id,nn,i
        integer ixdim,iydim,izdim,itdim
        character*80 name,time_v,time_d
        logical :: bdebug = .true.

        character(len=11), save :: xdims(4) =   (/
     +           'x          '
     +          ,'xpos       '
     +          ,'lon        '
     +          ,'west_east  '
     +                                          /)
        character(len=11), save :: ydims(4) =   (/
     +           'y          '
     +          ,'ypos       '
     +          ,'lat        '
     +          ,'south_north'
     +                                          /)
        character(len=15), save :: zdims(4) =   (/
     +           'z              '
     +          ,'zpos           '
     +          ,'bottom_top_stag'
     +          ,'level          '
     +                                          /)
        character(len=4), save :: tdims(2) =    (/
     +           'time'
     +          ,'Time'
     +                                          /)

        call nc_get_dim_totnum(ncid,ndims)

        ixdim = 0
        iydim = 0
        izdim = 0
        itdim = 0
        nx = 0
        ny = 0
        nz = 0
        nt = 0
        time_id = 0
        time_v = ' '

        do dim_id=1,ndims
          call nc_get_dim_name(ncid,dim_id,name)
          call nc_get_dim_len(ncid,dim_id,n)

          do i=1,size(xdims)
            if( name == xdims(i) ) ixdim = dim_id
          end do

          do i=1,size(ydims)
            if( name == ydims(i) ) iydim = dim_id
          end do

          do i=1,size(zdims)
            if( name == zdims(i) ) izdim = dim_id
          end do

          do i=1,size(tdims)
            if( name == tdims(i) ) itdim = dim_id
          end do
        end do

        if( bverb ) write(6,*) 'dimensions: '

        if( ixdim > 0 ) then
          call nc_get_dim_name(ncid,ixdim,name)
          call nc_get_dim_len(ncid,ixdim,nx)
          if( bverb ) write(6,*) '   xdim: ',nx,'  (',trim(name),')'
        end if

        if( iydim > 0 ) then
          call nc_get_dim_name(ncid,iydim,name)
          call nc_get_dim_len(ncid,iydim,ny)
          if( bverb ) write(6,*) '   ydim: ',ny,'  (',trim(name),')'
        end if

        if( izdim > 0 ) then
          call nc_get_dim_name(ncid,izdim,name)
          call nc_get_dim_len(ncid,izdim,nz)
          if( bverb ) write(6,*) '   zdim: ',nz,'  (',trim(name),')'
        end if

        if( itdim > 0 ) then    !handle time dimension
          call nc_get_dim_name(ncid,itdim,name)
          call nc_get_dim_len(ncid,itdim,nt)
          if( bverb ) write(6,*) '   tdim: ',nt,'  (',trim(name),')'
          time_d = name
          call nc_set_time_name(time_d,time_v)
        end if

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

        subroutine get_xycoord_names(ncid,bverb,xname,yname)

        implicit none

        integer ncid
        logical bverb
        character*(*) xname,yname

        integer nvars,iv,var_id
        character*80 name,atext

        xname = ' '
        yname = ' '

        call nc_get_var_totnum(ncid,nvars)

        do var_id=1,nvars

          call nc_get_var_name(ncid,var_id,name)

          call nc_get_var_attr(ncid,var_id,'standard_name',atext)
          if( atext == 'longitude' ) call set_name(xname,name)
          if( atext == 'latitude' ) call set_name(yname,name)

          call nc_get_var_attr(ncid,var_id,'long_name',atext)
          if( atext == 'longitude' ) call set_name(xname,name)
          if( atext == 'latitude' ) call set_name(yname,name)
          if( atext == 'Longitude' ) call set_name(xname,name)
          if( atext == 'Latitude' ) call set_name(yname,name)
	  if( atext == 'Longitude of scalars' ) call set_name(xname,name)
          if( atext == 'Latitude of scalars' ) call set_name(yname,name)

          call nc_get_var_attr(ncid,var_id,'description',atext)
          if( atext(1:10) == 'LONGITUDE,' ) call set_name(xname,name)
          if( atext(1:9) == 'LATITUDE,' ) call set_name(yname,name)

        end do

        if( bverb ) write(6,*) '   xcoord: ',trim(xname)
        if( bverb ) write(6,*) '   ycoord: ',trim(yname)

        end

c*****************************************************************

        subroutine get_tcoord_name(ncid,bverb,tcoord)

        implicit none

        integer ncid
        logical bverb
        character*(*) tcoord

        integer nvars,iv,var_id
        character*80 name,atext,time_d

        tcoord = ' '

        call nc_get_var_totnum(ncid,nvars)

        do var_id=1,nvars

          call nc_get_var_name(ncid,var_id,name)

          call nc_get_var_attr(ncid,var_id,'standard_name',atext)
          if( atext == 'time' ) call set_name(tcoord,name)

          call nc_get_var_attr(ncid,var_id,'long_name',atext)
          if( atext == 'time' ) call set_name(tcoord,name)
          if( atext == 'Julian day (UTC) of the station' ) 
     +				call set_name(tcoord,name)

          call nc_get_var_attr(ncid,var_id,'description',atext)
	  if( atext(1:13) == 'minutes since' ) call set_name(tcoord,name)

        end do

        time_d = ' '
        call nc_set_time_name(time_d,tcoord)

        if( bverb ) write(6,*) '   tcoord: ',trim(tcoord)

        end

c*****************************************************************

        subroutine get_zcoord_name(ncid,bverb,zcoord)

        implicit none

        integer ncid
        logical bverb
        character*(*) zcoord

        integer nvars,iv,var_id
        character*80 name,atext

        zcoord = ' '

        call nc_get_var_totnum(ncid,nvars)

        do var_id=1,nvars

          call nc_get_var_name(ncid,var_id,name)

          call nc_get_var_attr(ncid,var_id,'standard_name',atext)
          if( atext == 'zcoord' ) call set_name(zcoord,name)
          if( atext == 'sigma of cell face' ) call set_name(zcoord,name)

          call nc_get_var_attr(ncid,var_id,'long_name',atext)
          if( atext == 'zcoord' ) call set_name(zcoord,name)
          if( atext == 'sigma of cell face' ) call set_name(zcoord,name)

          call nc_get_var_attr(ncid,var_id,'description',atext)
          if( atext(1:18) == 'eta values on full' )
     +                          call set_name(zcoord,name)
          if( atext == 'bottom of vertical layers' )
     +                          call set_name(zcoord,name)
        end do

        if( bverb ) write(6,*) '   zcoord: ',trim(zcoord)

        end

c*****************************************************************

        subroutine set_name(varname,newname)

        implicit none

        character*(*) varname,newname

        if( varname == ' ' ) varname = newname

        end

c*****************************************************************
c*****************************************************************
c*****************************************************************

c       else if( iwhat .eq. 3 ) then    !wrf
c       'U10','u'  'V10','v'  'PSFC','p'
c       'RAINR','r',8
c       'T2','t'  'CLOUD','c'  'RH','h'  'SWDOWN','s'

c       else if( iwhat .eq. 4 ) then    !myocean
c       'sossheig','Z'  'vosaline','S'  'votemper','T'

c       else if( iwhat .eq. 5 ) then    !ROMS/TOMS
c       'salt','S' 'temp','T'

c       else if( iwhat .eq. 6 ) then    !wrf 2
c       'U10','u'  'V10','v'
c       'MSLP','p' 'PSFC','p'
c       'RAINR','r',rfact
c       'T2','t'  'CLOUD','c',0.01  'RH','h'  'SWDOWN','s'

c       else if( iwhat .eq. 7 ) then    !ECMWF
c       'var165','u'  'var166','v'  'var151','p'  'var228','r',rfact
c       'var167','t'  'var187','c'  'var168','d'  'var176','s'

c*****************************************************************
c*****************************************************************
c*****************************************************************


