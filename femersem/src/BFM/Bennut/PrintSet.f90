
SUBROUTINE PrintSet(NUTR)
     USE mem, ONLY:BoxNumberX,BoxNumberY
     USE global_mem, ONLY: LOGUNIT
     USE bennut_variables,ONLY:sets
     implicit none
     integer  ::NUTR    !Specification`

     integer                      ::i,k,l,n
    
     write(LOGUNIT,'(''Warnings benthic nutrient model'')')
     write(LOGUNIT,'(''Point x,y:'',I5,'','',I5 )') BoxNumberX,BoxNumberY
     write(LOGUNIT,'('' Layer Definition'')')

     do i= 1,sets(NUTR)%equa
       write( LOGUNIT,'(''layer'',I6,G12.5,G12.5,G12.5)') i,  &
              sets(NUTR)%diff(i), sets(NUTR)%poro(i), sets(NUTR)%ads(i)
     enddo
     write(LOGUNIT,*)

     write(LOGUNIT,'(''"Set Definition'')')
     do i= 1,sets(NUTR)%nn
           k=sets(NUTR)%coeffs(i)%il
           l=sets(NUTR)%coeffs(i)%ia
           if ( l>=0)  then
              write(LOGUNIT,'(I5,'' term:'',I2,'' type:'',i2)') i, k,l
           else 
              write(LOGUNIT, &
                '(I5,'' term:'',I2,'' type:'',i2,'' lambda:'',2G12.5)') i,k,l, &
                 sets(NUTR)%coeffs%labda(1),sets(NUTR)%coeffs%labda(2)
           endif
     enddo
     write(LOGUNIT,*)
END 
