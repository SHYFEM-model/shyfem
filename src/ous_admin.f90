!
! $Id: subousa.f,v 1.2 2005/11/03 16:59:25 georg Exp $
!
! OUS file administration routines
!
! contents :
!
! subroutine wrousa
!
! revision log :
!
! 26.01.1998	ggu	$$ITMOUT - adjust itmout for first write
! 01.09.2003	ggu	new routine wrousa
! 02.09.2003	ggu	bug fix in wrousa: save nbout
! 25.11.2004	ggu	in new file subousa.f
! 18.05.2005	ggu	initial itmout is changed (now itanf)
! 20.01.2014	ggu	new calls for ous writing implemented
! 15.10.2015	ggu	added new calls for shy file format
!
!********************************************************
!-----------------------------------------------------------------
        module ous_admin
!-----------------------------------------------------------------
        contains
!-----------------------------------------------------------------

	subroutine wrousa

! writes and administers ous file

        use depth
        use hydro_admin
        use levels
        use shympi
        use basin, only : nkn,nel,ngr,mbw,nkndi,neldi
        use mpi_io_admin
        use defnames
        use para
        use output
        use version
        use outputd
        use ioous
        use shy_util

        implicit none

        include 'femtime.h'
        include 'simul.h'

        integer itmout,ierr
        double precision href,hzoff
        integer date,time
        character*80 title,femver

        integer iround
        integer wfout,wrout
        integer nvar,ftype,id
        double precision dtime
        character*80 file

        integer idtout,itout
        integer icall,nbout,nvers
        integer ishyff
        integer, save :: ia_out(4)
        double precision, save :: da_out(4)
        save idtout,itout
        save icall,nvers,nbout
        data icall,nvers,nbout /0,2,0/

        ishyff = nint(getpar('ishyff'))

        if( icall .eq. -1 ) return

        if(shympi_partition_on_nodes()) return

        if(bmpi) then
          if( icall .eq. 0 ) then
            call rebuild_ous_header    
          end if
        end if

        if( icall .eq. 0 ) then
                ia_out = 0
                da_out = 0.
                call init_output('itmout','idtout',ia_out)
                if( ishyff == 1 ) ia_out = 0
                call init_output_d('itmout','idtout',da_out)
                if( ishyff == 0 ) da_out = 0

                if( .not. has_output(ia_out) .and. .not. has_output_d(da_out) ) icall = -1
                if( icall .eq. -1 ) return
                
		!if( shympi_is_parallel() ) then
		!  stop 'error stop wrousa: not mpi ready'
		!end if

                if( has_output_d(da_out) ) then
                  nvar = 4
                  ftype = 1
                  call shy_make_output_name('.hydro.shy',file)
                  call shy_open_output_file(file,3,nlv,nvar,ftype,id)
                  da_out(4) = id
                end if

                if( has_output(ia_out) ) then

                nbout = ifemop('.ous','unformatted','new')
                if(nbout.le.0) goto 77
                ia_out(4) = nbout

                if(shympi_is_master()) then

                href=getpar('href')             !reference level
                hzoff=getpar('hzoff')           !minimum depth
                date = nint(dgetpar('date'))
                time = nint(dgetpar('time'))
                title = descrp
                call get_shyfem_version(femver)

                call ous_init(nbout,nvers)
                call ous_set_title(nbout,title)
                call ous_set_date(nbout,date,time)
                call ous_set_femver(nbout,femver)
                call ous_set_hparams(nbout,href,hzoff)
                call ous_write_header(nbout,nkndi,neldi,nlvdi,ierr)
                if(ierr.gt.0) goto 78
                if(bmpi) then
                  call ous_write_header2(nbout,outIlhv,hlv,outHev,ierr) 
                else
                  call ous_write_header2(nbout,ilhv,hlv,hev,ierr)
                end if
                if(ierr.gt.0) goto 75

                end if

                end if
        end if

        icall = icall + 1

        if( next_output(ia_out) ) then
          if(bmpi) then
            call rebuild_structures(ilhv,znv,zenv,utlnv,vtlnv)
            if(shympi_is_master()) then
             call ous_write_record(nbout,it,nlvdi,outIlhv,outZnv,outZenv,outUtlnv,outVtlnv,ierr)
             if(ierr.ne.0.) goto 79
            end if
          else
          call ous_write_record(nbout,it,nlvdi,ilhv,znv,zenv,utlnv,vtlnv,ierr)
          if(ierr.ne.0.) goto 79
          end if
        end if

        if( next_output_d(da_out) ) then
          if(bmpi) then
            call rebuild_structures(ilhv,znv,zenv,utlnv,vtlnv)
          end if
          if(shympi_is_master()) then
            id = nint(da_out(4))
            dtime = t_act
            call shy_write_hydro_records(id,dtime,nlvdi,outZnv,outZenv,outUtlnv,outVtlnv)
          end if
        end if

        return
   77   continue
        write(6,*) 'Error opening OUS file :'
        stop 'error stop : wrousa'
   78   continue
        write(6,*) 'Error writing first header of OUS file'
        write(6,*) 'unit,err :',nbout,ierr
        stop 'error stop : wrousa'
   75   continue
        write(6,*) 'Error writing second header of OUS file'
        write(6,*) 'unit,err :',nbout,ierr
        stop 'error stop : wrousa'
   79   continue
        write(6,*) 'Error writing data record of OUS file'
        write(6,*) 'unit,err :',nbout,ierr
        stop 'error stop : wrousa'
        end

!********************************************************

!-----------------------------------------------------------------
        end module ous_admin
!-----------------------------------------------------------------
