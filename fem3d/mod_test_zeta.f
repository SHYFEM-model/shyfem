
!==================================================================
        module mod_test_zeta
!==================================================================

        implicit none

	integer, private, save :: nn_step = 0

!==================================================================
        contains
!==================================================================

	subroutine test_zeta_init
	use basin
	use shympi

	implicit none

        character*4 my_id_s,n_threads_s
	real azpar,ampar

#if !defined(test_zeta)
	return
#endif

        write(n_threads_s,'(i4.4)')n_threads
        write(my_id_s,'(i4.4)')my_id+1

	call getazam(azpar,ampar)

	if( (azpar == 0. .and. ampar == 1. ) .or.
     +        (azpar == 1. .and. ampar == 0. ) )  then
           open(unit = 99998, 
     +          file = "test_zeta."//trim(n_threads_s)//
     +                 ".EXPLICIT_"//trim(my_id_s))
	else 
           open(unit = 99998, 
     +          file = "test_zeta."//trim(n_threads_s)//
     +                 ".IMPLICIT_"//trim(my_id_s))
        endif

        if(nkn<100) then
          nn_step=1
        else
          nn_step=nkn/100
        endif

	end subroutine test_zeta_init

!******************************************************************

	subroutine test_zeta_write

	use basin
	use mod_hydro

	implicit none

	integer nn

#if !defined(test_zeta)
	return
#endif

        write(99998,'(i10)') 0   
        do nn=1,nkn,nn_step 
          write(99998,'(i10,E25.15)') nn,znv(nn)   
        end do 

	end subroutine test_zeta_write

!******************************************************************

	subroutine test_zeta_close

	implicit none

#if !defined(test_zeta)
	return
#endif

	close(99998)

	end subroutine test_zeta_close

!==================================================================
        end module mod_test_zeta
!==================================================================

