!************************************************************************
!**************************co2sys routines***************************
!************************************************************************
!Includes:
!  module AQUABC_II_GLOBAL
!  module VECTOR_MATRIX_UTILS
!  module CO2SYS_CDIAC


 module AQUABC_II_GLOBAL
    implicit none
    integer :: DBL_PREC
    parameter(DBL_PREC = selected_real_kind(15, 307))
 end module AQUABC_II_GLOBAL

!*****************************************************
!*****************************************************

 module VECTOR_MATRIX_UTILS
    use AQUABC_II_GLOBAL
    implicit none
 contains
    !**************************************************************************
    ! UTILITY FUNCTIONS BETWEEN FORTRAN AND MATLAB
    !**************************************************************************
    function SAME_ELEMENTS_INETEGER_VECTOR(VECTOR) result(UNIQUE)

        !INGOING VARIABLES
        real(kind = DBL_PREC), intent(in), dimension(:) :: VECTOR
        !END OF INGOING VARIABLES

        !RETURNED VARIABLE
        integer :: UNIQUE
        !END OF RETURNED VARIABLE

        ! AUXILLARY VARIABLES
        real(kind = DBL_PREC) :: FIRST_ELEMENTS_VALUE
        integer :: i
        ! END OF AUXILLARY VARIABLES

        UNIQUE = 1
        FIRST_ELEMENTS_VALUE = VECTOR(1)

        do i = 1, size(VECTOR)
            if (VECTOR(i).ne.FIRST_ELEMENTS_VALUE) then
                UNIQUE = 0
                return
            endif
        enddo

    end function SAME_ELEMENTS_INETEGER_VECTOR



    subroutine ASSIGN_DBL_VECTOR_CONTENT(ARRAY_TO, ARRAY_FROM)
    ! subroutine is not used more?
    
        !INGOING VARIABLES
        real(kind = DBL_PREC), intent(in)   , dimension(:) :: ARRAY_FROM
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: ARRAY_TO
        !END OF INGOING VARIABLES

        if(allocated(ARRAY_TO)) then
            deallocate(ARRAY_TO)
        endif

        allocate(ARRAY_TO(size(ARRAY_FROM)))
        ARRAY_TO(:) = ARRAY_FROM(:)
    end subroutine ASSIGN_DBL_VECTOR_CONTENT


    subroutine ASSIGN_INT_VECTOR_CONTENT(ARRAY_TO, ARRAY_FROM)
     ! subroutine is not used more?
        !INGOING VARIABLES
        integer, intent(in)   , dimension(:) :: ARRAY_FROM
        integer, intent(inout), allocatable, dimension(:) :: ARRAY_TO
        !END OF INGOING VARIABLES

        if(allocated(ARRAY_TO)) then
            deallocate(ARRAY_TO)
        endif

        allocate(ARRAY_TO(size(ARRAY_FROM)))
        ARRAY_TO(:) = ARRAY_FROM(:)
    end subroutine ASSIGN_INT_VECTOR_CONTENT


    subroutine LOGICAL_VECTOR_TO_INT_VECTOR(LOGICAL_VECTOR, INT_VECTOR)
    ! subroutine is not used more?
        logical, intent(in)   , dimension(:)              :: LOGICAL_VECTOR
        integer, intent(inout), allocatable, dimension(:) :: INT_VECTOR
        integer :: i

        if(size(LOGICAL_VECTOR) > 0) then
            if(allocated(INT_VECTOR)) then
                deallocate(INT_VECTOR)
            end if

        allocate(INT_VECTOR(size(LOGICAL_VECTOR)))
        
            INT_VECTOR = 0

            where(LOGICAL_VECTOR)
                INT_VECTOR = 1
            end where
        else
            write(*,*) 'Error in subroutine '
            write(*,*) 'LOGICAL_VECTOR_TO_INT_VECTOR(LOGICAL_VECTOR, INT_VECTOR)'
            write(*,*) 'The size of the argument LOGICAL_VECTOR passed to subroutine '
            write(*,*) 'LOGICAL_VECTOR_TO_INT_VECTOR should be greater than zero'
            stop            
        end if
    end subroutine LOGICAL_VECTOR_TO_INT_VECTOR


    subroutine GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)
        integer, intent(in)   , allocatable, dimension(:) :: INDEX_ARRAY_MASK
        integer, intent(inout), allocatable, dimension(:) :: INDEX_ARRAY
        integer :: i
        integer :: j

        if (allocated(INDEX_ARRAY)) then
            deallocate(INDEX_ARRAY)
        end if

        allocate(INDEX_ARRAY(sum(INDEX_ARRAY_MASK)))
        j = 1

        do i = lbound(INDEX_ARRAY_MASK,1), ubound(INDEX_ARRAY_MASK,1)
            if(INDEX_ARRAY_MASK(i) == 1) then
                INDEX_ARRAY(j) = i
                j = j + 1
            end if
        end do

    end subroutine GENERATE_INDEX_ARRAY

    ! END OF UTILITY FUNCTIONS BETWEEN FORTAN AND MATLAB
 end module VECTOR_MATRIX_UTILS

!*****************************************************
!*****************************************************




!*****************************************************
!*****************************************************

 module CO2SYS_CDIAC
    use AQUABC_II_GLOBAL
    use VECTOR_MATRIX_UTILS
    implicit none

 contains

    subroutine CO2SYS &
               (PAR1     , PAR2  , PAR1TYPE, PAR2TYPE, SALT, TEMPIN,   &
                TEMPOUT  , PRESIN, PRESOUT , SI , PO4, pHSCALEIN  ,   &
                K1K2CONSTANTS, KSO4CONSTANTS, &
                CO2SYS_OUT_DATA, NICEHEADERS, &
                ntps)

        !INGOING VARIABLES
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: PAR1
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: PAR2
!         integer              , intent(in)   , allocatable, dimension(:)   :: PAR1TYPE
!         integer              , intent(in)   , allocatable, dimension(:)   :: PAR2TYPE
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: SALT
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: TEMPIN
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: TEMPOUT
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: PRESIN
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: PRESOUT
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: SI
!         real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:)   :: PO4
!         integer              , intent(in)   , allocatable, dimension(:)   :: pHSCALEIN
!         integer              , intent(in)   , allocatable, dimension(:)   :: K1K2CONSTANTS
!         integer              , intent(in)   , allocatable, dimension(:)   :: KSO4CONSTANTS
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: PAR1
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: PAR2
         integer              , intent(in),    dimension(ntps)   :: PAR1TYPE
         integer              , intent(in),    dimension(ntps)   :: PAR2TYPE
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: SALT
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: TEMPIN
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: TEMPOUT
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: PRESIN
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: PRESOUT
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: SI
         real(kind = DBL_PREC), intent(in),    dimension(ntps)   :: PO4
         integer              , intent(in),    dimension(ntps)   :: pHSCALEIN
         integer              , intent(in),    dimension(ntps)   :: K1K2CONSTANTS
         integer              , intent(in),    dimension(ntps)   :: KSO4CONSTANTS         
        
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:,:) :: CO2SYS_OUT_DATA
        character(len = 34)  , intent(inout), allocatable, dimension(:)   :: NICEHEADERS
        !END OF INGOING VARIABLES

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, allocatable, dimension(:) :: pHScale
        integer, allocatable, dimension(:) :: WhichKs
        integer, allocatable, dimension(:) :: WhoseKSO4
        integer, allocatable, dimension(:) :: p1
        integer, allocatable, dimension(:) :: p2

        real(kind = DBL_PREC), allocatable, dimension(:) :: TempCi
        real(kind = DBL_PREC), allocatable, dimension(:) :: TempCo
        real(kind = DBL_PREC), allocatable, dimension(:) :: Pdbari
        real(kind = DBL_PREC), allocatable, dimension(:) :: Pdbaro
        real(kind = DBL_PREC), allocatable, dimension(:) :: Sal
        real(kind = DBL_PREC), allocatable, dimension(:) :: sqrSal
        real(kind = DBL_PREC), allocatable, dimension(:) :: TP
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), allocatable, dimension(:) :: TB
        real(kind = DBL_PREC), allocatable, dimension(:) :: TF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TS

        real(kind = DBL_PREC), allocatable, dimension(:) :: TempK
        real(kind = DBL_PREC), allocatable, dimension(:) :: RT
        real(kind = DBL_PREC), allocatable, dimension(:) :: logTempK
        real(kind = DBL_PREC), allocatable, dimension(:) :: Pbar
        real(kind = DBL_PREC), allocatable, dimension(:) :: TempK100
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnK0
        real(kind = DBL_PREC), allocatable, dimension(:) :: K0
        real(kind = DBL_PREC), allocatable, dimension(:) :: KS

        real(kind = DBL_PREC), allocatable, dimension(:) :: fH
        real(kind = DBL_PREC), allocatable, dimension(:) :: K1
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2
        real(kind = DBL_PREC), allocatable, dimension(:) :: KW
        real(kind = DBL_PREC), allocatable, dimension(:) :: KB
        real(kind = DBL_PREC), allocatable, dimension(:) :: KF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSi

        real(kind = DBL_PREC), allocatable, dimension(:) :: FugFac
        real(kind = DBL_PREC), allocatable, dimension(:) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        integer :: ERROR

        ! AUXILLARY VARIABLES
        real(kind = DBL_PREC), dimension(1:14) :: VECLENGTHS

        real(kind = DBL_PREC), allocatable, dimension(:) :: TA    ! Talk
        real(kind = DBL_PREC), allocatable, dimension(:) :: TC    ! DIC
        real(kind = DBL_PREC), allocatable, dimension(:) :: PH    ! pH
        real(kind = DBL_PREC), allocatable, dimension(:) :: PC    ! pCO2
        real(kind = DBL_PREC), allocatable, dimension(:) :: FC    ! fCO2
        real(kind = DBL_PREC), allocatable, dimension(:) :: PengCorrection
        real(kind = DBL_PREC), allocatable, dimension(:) :: TAc
        real(kind = DBL_PREC), allocatable, dimension(:) :: TAc_Interface
        real(kind = DBL_PREC), allocatable, dimension(:) :: TAc_minus_Peng
        real(kind = DBL_PREC), allocatable, dimension(:) :: TCc
        real(kind = DBL_PREC), allocatable, dimension(:) :: TCc_Interface
        real(kind = DBL_PREC), allocatable, dimension(:) :: PHic
        real(kind = DBL_PREC), allocatable, dimension(:) :: PHic_Interface
        real(kind = DBL_PREC), allocatable, dimension(:) :: PCic
        !real(kind = DBL_PREC), allocatable, dimension(:) :: PCic_Interface
        real(kind = DBL_PREC), allocatable, dimension(:) :: FCic
        real(kind = DBL_PREC), allocatable, dimension(:) :: FCic_Interface

        integer, allocatable, dimension(:) :: Icase
        integer, allocatable, dimension(:) :: INDEX_ARRAY
        integer, allocatable, dimension(:) :: INDEX_ARRAY_MASK

        real(kind = DBL_PREC), allocatable, dimension(:) :: CO2inp
        real(kind = DBL_PREC), allocatable, dimension(:) :: HCO3inp
        real(kind = DBL_PREC), allocatable, dimension(:) :: CO3inp
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlkinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: OHinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlkinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlkinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfreeinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4inp
        real(kind = DBL_PREC), allocatable, dimension(:) :: HFinp

        real(kind = DBL_PREC), allocatable, dimension(:) :: Revelleinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: OmegaCainp
        real(kind = DBL_PREC), allocatable, dimension(:) :: OmegaArinp
        real(kind = DBL_PREC), allocatable, dimension(:) :: xCO2dryinp

        real(kind = DBL_PREC), allocatable, dimension(:) :: pHicT
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHicS
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHicF
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHicN

        real(kind = DBL_PREC), allocatable, dimension(:,:) :: KIVEC
        real(kind = DBL_PREC), allocatable, dimension(:,:) :: KOVEC
        real(kind = DBL_PREC), allocatable, dimension(:,:) :: TVEC

        real(kind = DBL_PREC), allocatable, dimension(:) :: pHoc
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHocT
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHocS
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHocF
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHocN

        real(kind = DBL_PREC), allocatable, dimension(:) :: FCoc
        real(kind = DBL_PREC), allocatable, dimension(:) :: PCoc

        real(kind = DBL_PREC), allocatable, dimension(:) :: CO2out
        real(kind = DBL_PREC), allocatable, dimension(:) :: HCO3out
        real(kind = DBL_PREC), allocatable, dimension(:) :: CO3out
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlkout
        real(kind = DBL_PREC), allocatable, dimension(:) :: OHout
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlkout
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlkout
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfreeout
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4out
        real(kind = DBL_PREC), allocatable, dimension(:) :: HFout

        real(kind = DBL_PREC), allocatable, dimension(:) :: Revelleout
        real(kind = DBL_PREC), allocatable, dimension(:) :: OmegaCaout
        real(kind = DBL_PREC), allocatable, dimension(:) :: OmegaArout
        real(kind = DBL_PREC), allocatable, dimension(:) :: xCO2dryout
        ! END OF AUXILLARY VARIABLES

        integer :: DEBUG
        parameter(DEBUG = 0)
        


        ERROR = 0

        VECLENGTHS(1)  = size(PAR1)
        VECLENGTHS(2)  = size(PAR2)
        VECLENGTHS(3)  = size(PAR1TYPE)
        VECLENGTHS(4)  = size(PAR2TYPE)
        VECLENGTHS(5)  = size(SALT)
        VECLENGTHS(6)  = size(TEMPIN)
        VECLENGTHS(7)  = size(TEMPOUT)
        VECLENGTHS(8)  = size(PRESIN)
        VECLENGTHS(9)  = size(PRESOUT)
        VECLENGTHS(10) = size(SI)
        VECLENGTHS(11) = size(PO4)
        VECLENGTHS(12) = size(pHSCALEIN)
        VECLENGTHS(13) = size(K1K2CONSTANTS)
        VECLENGTHS(14) = size(KSO4CONSTANTS)

        ntps = maxval(VECLENGTHS)

        if (SAME_ELEMENTS_INETEGER_VECTOR(VECLENGTHS) == 0) then
            write(unit = *, fmt = *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit = *, fmt = *) '!!!    ERROR in function CO2SYS    !!!'
            write(unit = *, fmt = *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(unit = *, fmt = *) ''
            write(unit = *, fmt = *) 'All the input arguments must have the same number of elements'
            write(unit = *, fmt = *) ''
            write(unit = *, fmt = '(a47, i10)') 'Argument : PAR1         , Number of elements = ', VECLENGTHS(1)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PAR2         , Number of elements = ', VECLENGTHS(2)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PAR1TYPE     , Number of elements = ', VECLENGTHS(3)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PAR2TYPE     , Number of elements = ', VECLENGTHS(4)
            write(unit = *, fmt = '(a47, i10)') 'Argument : SAL          , Number of elements = ', VECLENGTHS(5)
            write(unit = *, fmt = '(a47, i10)') 'Argument : TEMPIN       , Number of elements = ', VECLENGTHS(6)
            write(unit = *, fmt = '(a47, i10)') 'Argument : TEMPOUT      , Number of elements = ', VECLENGTHS(7)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PRESIN       , Number of elements = ', VECLENGTHS(8)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PRESOUT      , Number of elements = ', VECLENGTHS(9)
            write(unit = *, fmt = '(a47, i10)') 'Argument : SI           , Number of elements = ', VECLENGTHS(10)
            write(unit = *, fmt = '(a47, i10)') 'Argument : PO4          , Number of elements = ', VECLENGTHS(11)
            write(unit = *, fmt = '(a47, i10)') 'Argument : pHSCALEIN    , Number of elements = ', VECLENGTHS(12)
            write(unit = *, fmt = '(a47, i10)') 'Argument : K1K2CONSTANTS, Number of elements = ', VECLENGTHS(13)
            write(unit = *, fmt = '(a47, i10)') 'Argument : KSO4CONSTANTS, Number of elements = ', VECLENGTHS(14)
            stop
        endif

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #1'
        endif

        !Allocate all the former CO2SYS global variables
        allocate(pHScale  (ntps))
        allocate(WhichKs  (ntps))
        allocate(WhoseKSO4(ntps))
        allocate(p1       (ntps))
        allocate(p2       (ntps))
        allocate(TempCi   (ntps))
        allocate(TempCo   (ntps))
        allocate(Pdbari   (ntps))
        allocate(Pdbaro   (ntps))
        allocate(Sal      (ntps))
        allocate(sqrSal   (ntps))
        allocate(TP       (ntps))
        allocate(TSi      (ntps))
        allocate(TB       (ntps))
        allocate(TF       (ntps))
        allocate(TS       (ntps))
        allocate(TempK    (ntps))
        allocate(RT       (ntps))
        allocate(logTempK (ntps))
        allocate(Pbar     (ntps))
        allocate(TempK100 (ntps))
        allocate(lnK0     (ntps))
        allocate(K0       (ntps))
        allocate(KS       (ntps))
        allocate(fH       (ntps))
        allocate(K1       (ntps))
        allocate(K2       (ntps))
        allocate(KW       (ntps))
        allocate(KB       (ntps))
        allocate(KF       (ntps))
        allocate(KP1      (ntps))
        allocate(KP2      (ntps))
        allocate(KP3      (ntps))
        allocate(KSi      (ntps))
        allocate(FugFac   (ntps))
        allocate(VPFac    (ntps))
        !End of allocate all the former CO2SYS global variables
        
        
        !call ASSIGN_INT_VECTOR_CONTENT(pHScale  , pHSCALEIN)     ! MATLAB CODE :  pHScale      = pHSCALEIN;
        if(.not. allocated(pHScale))      allocate(pHScale(ntps))
        pHScale  = pHSCALEIN
        !call ASSIGN_INT_VECTOR_CONTENT(WhichKs  , K1K2CONSTANTS) ! MATLAB CODE :  WhichKs      = K1K2CONSTANTS;
        if(.not. allocated(WhichKs))      allocate(WhichKs(ntps))        
        WhichKs  = K1K2CONSTANTS
        !call ASSIGN_INT_VECTOR_CONTENT(WhoseKSO4, KSO4CONSTANTS) ! MATLAB CODE :  WhoseKSO4    = KSO4CONSTANTS;
        if(.not. allocated(WhoseKSO4))      allocate(WhoseKSO4(ntps))        
        WhoseKSO4 = KSO4CONSTANTS
        !call ASSIGN_INT_VECTOR_CONTENT(p1       , PAR1TYPE)      ! MATLAB CODE :  p1           = PAR1TYPE;
        if(.not. allocated(p1))      allocate(p1(ntps))        
        p1  = PAR1TYPE
        !call ASSIGN_INT_VECTOR_CONTENT(p2       , PAR2TYPE)      ! MATLAB CODE :  p2           = PAR2TYPE;
        if(.not. allocated(p2))      allocate(p2(ntps ))        
        p2 = PAR2TYPE
        !call ASSIGN_DBL_VECTOR_CONTENT(TempCi   , TEMPIN)        ! MATLAB CODE :  TempCi       = TEMPIN;
        if(.not. allocated(TempCi))      allocate(TempCi(ntps))        
        TempCi  = TEMPIN
        !call ASSIGN_DBL_VECTOR_CONTENT(TempCo   , TEMPOUT)       ! MATLAB CODE :  TempCo       = TEMPOUT;
        if(.not. allocated(TempCo))      allocate(TempCo(ntps))        
        TempCo  = TEMPOUT
        !call ASSIGN_DBL_VECTOR_CONTENT(Pdbari   , PRESIN)        ! MATLAB CODE :  Pdbari       = PRESIN;
        if(.not. allocated(Pdbari))      allocate(Pdbari(ntps))        
        Pdbari = PRESIN
        !call ASSIGN_DBL_VECTOR_CONTENT(Pdbaro   , PRESOUT)       ! MATLAB CODE :  Pdbaro       = PRESOUT;
        if(.not. allocated(Pdbaro))      allocate(Pdbaro(ntps))        
        Pdbaro  = PRESOUT
        !call ASSIGN_DBL_VECTOR_CONTENT(Sal      , SALT)          ! MATLAB CODE :  Sal          = SAL;
        if(.not. allocated(Sal))      allocate(Sal(ntps))        
        Sal  = SALT
        !call ASSIGN_DBL_VECTOR_CONTENT(sqrSal   , sqrt(SALT))    ! MATLAB CODE :  sqrSal       = sqrt(SAL);
        if(.not. allocated(sqrSal))      allocate(sqrSal(ntps))        
        sqrSal = sqrt(SALT)
        !call ASSIGN_DBL_VECTOR_CONTENT(TP       , PO4)           ! MATLAB CODE :  TP           = PO4;
        if(.not. allocated(TP ))      allocate(TP(ntps))        
        TP  = PO4
        !call ASSIGN_DBL_VECTOR_CONTENT(TSi      , SI)            ! MATLAB CODE :  TSi          = SI;
        if(.not. allocated(TSi ))      allocate(TSi(ntps))        
        TSi = SI


        
        RGasConstant = 83.1451D0  !RGasConstant = 83.1451;  % ml bar-1 K-1 mol-1, DOEv2

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #2'
        endif

        ! Generate empty vectors for...
        allocate(TA(ntps), &    ! Talk
                 TC(ntps), &    ! DIC
                 PH(ntps), &    ! pH
                 PC(ntps), &    ! pCO2
                 FC(ntps), &     ! fCO2
                 PengCorrection(ntps))

        TA             = 0.0D0
        TC             = 0.0D0
        PH             = 0.0D0
        PC             = 0.0D0
        FC             = 0.0D0
        PengCorrection = 0.0D0

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #3'
        endif

        ! Assign values to empty vectors.

        ! Convert from micromol/kg to mol/kg
        ! MATLAB CODE :  F=(p1==1); TA(F)=PAR1(F)/1e6;
        where (p1==1)
            TA = PAR1 / 1.0D6
        end where

        ! Convert from micromol/kg to mol/kg
        ! MATLAB CODE :  F=(p1==2); TC(F)=PAR1(F)/1e6;
        where (p1==2)
            TC = PAR1 / 1.0D6
        end where

        ! MATLAB CODE :  F=(p1==3); PH(F)=PAR1(F);
        where (p1==3)
            PH = PAR1
        end where

        ! Convert from microatm. to atm.
        ! MATLAB CODE :  F=(p1==4); PC(F)=PAR1(F)/1e6;
        where (p1==4)
            PC = PAR1 / 1.0D6
        end where

        ! Convert from microatm. to atm.
        ! MATLAB CODE :  F=(p1==5); FC(F)=PAR1(F)/1e6;
        where (p1==5)
            FC = PAR1 / 1.0D6
        end where

        ! Convert from micromol/kg to mol/kg
        ! MATLAB CODE :   F=(p2==1); TA(F)=PAR2(F)/1e6;
        where (p2==1)
            TA = PAR2 / 1.0D6
        end where

        ! Convert from micromol/kg to mol/kg
        ! MATLAB CODE :   F=(p2==2); TC(F)=PAR2(F)/1e6;
        where (p2==2)
            TC = PAR2 / 1.0D6
        end where

        ! MATLAB CODE :  F=(p1==3); PH(F)=PAR1(F);
        where (p1==3)
            PH = PAR1
        end where

        ! Convert from microatm. to atm.
        ! MATLAB CODE :  F=(p2==4); PC(F)=PAR2(F)/1e6;
        where (p1==4)
            PC = PAR2 / 1.0D6
        end where

        ! Convert from microatm. to atm.
        ! MATLAB CODE :  F=(p2==5); FC(F)=PAR2(F)/1e6;
        where (p2==5)
            FC = PAR2 / 1.0D6
        end where

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #4'
        endif

        ! All other cases
        ! The following code was origially after conditions
        ! F=(WhichKs==8); and F=(WhichKs==8 | WhichKs==6);
        ! MATLAB CODE :  F=~F;
        ! MATLAB CODE :  TP(F)  = TP(F)./1e6;
        ! MATLAB CODE :  TSi(F) = TSi(F)./1e6;
        TP  = TP  / 1.0D6
        TSi = TSi / 1.0D6

        ! Generate the columns holding Si, Phos and Sal.
        ! Pure Water case:
        ! MATLAB CODE :  F=(WhichKs==8);
        ! MATLAB CODE :  Sal(F) = 0;
        where (WhichKs==8)
            Sal = 0.0D0
        end where

        ! GEOSECS and Pure Water:
        ! MATLAB CODE :  F=(WhichKs==8 | WhichKs==6);
        ! MATLAB CODE :  TP(F)  = 0;
        ! MATLAB CODE :  TSi(F) = 0;
        where (WhichKs==8 .or. WhichKs==6)
            TP  = 0.0D0
            TSi = 0.0D0
        end where

        ! The vector 'PengCorrection' is used to modify the value of TA, for those
        ! cases where WhichKs==7, since PAlk(Peng) = PAlk(Dickson) + TP.
        ! Thus, PengCorrection is 0 for all cases where WhichKs is not 7
        ! MATLAB CODE :  PengCorrection=zeros(ntps,1); F=WhichKs==7; PengCorrection(F)=TP(F);
        where (WhichKs==7)
            PengCorrection = TP
        end where

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #5'
        endif

        ! Calculate the constants for all samples at input conditions
        ! The constants calculated for each sample will be on the appropriate pH scale!
        ! MATLAB CODE :  Constants(TempCi,Pdbari);
        

        
        call Constants &
             (TempCi,Pdbari, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)
              


        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #6'
        endif

        ! Make sure fCO2 is available for each sample that has pCO2.
        !MATLAB CODE: F=find(p1==4 | p2==4); FC(F) = PC(F).*FugFac(F);
        where((p1==4) .or. (p2==4))
            FC = PC * FugFac
        end where

        ! Generate vector for results, and copy the raw input values into them. This
        ! copies ~60% NaNs, which will be replaced for calculated values later on.
        call ASSIGN_DBL_VECTOR_CONTENT(TAc , TA)       !MATLAB CODE: TAc  = TA;
        call ASSIGN_DBL_VECTOR_CONTENT(TCc , TC)       !MATLAB CODE: TCc  = TC;
        call ASSIGN_DBL_VECTOR_CONTENT(PHic, PH)       !MATLAB CODE: PHic = PH;
        call ASSIGN_DBL_VECTOR_CONTENT(PCic, PC)       !MATLAB CODE: PCic = PC;
        call ASSIGN_DBL_VECTOR_CONTENT(FCic, FC)       !MATLAB CODE: FCic = FC;

        call ASSIGN_INT_VECTOR_CONTENT(Icase, 10 * min(p1,p2) + max(p1,p2))

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7'
        endif

        ! Calculate missing values for AT,CT,PH,FC:
        ! pCO2 will be calculated later on, routines work with fCO2.

        allocate(INDEX_ARRAY_MASK(ntps))

        ! *******************************************************************************
        ! input TA, TC
        ! *******************************************************************************

        ! Changes were needed to adapt the following MATLAB code to fortran
        ! F=Icase==12; % input TA, TC
        ! if any(F)
             ![PHic(F) FCic(F)] = CalculatepHfCO2fromTATC(TAc(F)-PengCorrection(F), TCc(F));
        !end

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-1'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==12)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            allocate(TAc_minus_Peng(sum(INDEX_ARRAY_MASK)), TCc_Interface (sum(INDEX_ARRAY_MASK)), &
                     PHic_Interface(sum(INDEX_ARRAY_MASK)), FCic_Interface(sum(INDEX_ARRAY_MASK)))

            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            TAc_minus_Peng = TAc(INDEX_ARRAY) - PengCorrection(INDEX_ARRAY)
            TCc_Interface  = TCc(INDEX_ARRAY)

            call CalculatepHfCO2fromTATC &
                 (DEBUG, TAc_minus_Peng , TCc_Interface, PHic_Interface, FCic_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            PHic(INDEX_ARRAY) = PHic_Interface
            FCic(INDEX_ARRAY) = FCic_Interface

            deallocate(TAc_minus_Peng, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of input TA, TC
        ! *******************************************************************************


        ! *******************************************************************************
        ! input TA, pH
        ! *******************************************************************************
        ! Changes were needed to adapt the following MATLAB code to fortran
        !F=Icase==13; % input TA, pH
        !if any(F)
        !    TCc(F)  =   CalculateTCfromTApH(TAc(F)-PengCorrection(F), PHic(F));
        !    FCic(F) = CalculatefCO2fromTCpH(TCc(F), PHic(F));
        !end
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-2'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==13)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            allocate(TAc_minus_Peng(sum(INDEX_ARRAY)), TCc_Interface (sum(INDEX_ARRAY)), &
                     PHic_Interface(sum(INDEX_ARRAY)), FCic_Interface(sum(INDEX_ARRAY)))

            TAc_minus_Peng = TAc(INDEX_ARRAY) - PengCorrection(INDEX_ARRAY)
            PHic_Interface = PHic(INDEX_ARRAY)

            call CalculateTCfromTApH  &
                 (TAc_minus_Peng, PHic_Interface, TCc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            call CalculatefCO2fromTCpH &
                 (TCc_Interface , PHic_Interface, FCic_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            TCc (INDEX_ARRAY) = TCc_Interface
            FCic(INDEX_ARRAY) = FCic_Interface

            deallocate(TAc_minus_Peng, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of  input TA, pH
        ! *******************************************************************************


        ! *******************************************************************************
        ! input TA, (pCO2 or fCO2)
        ! *******************************************************************************
        ! Changes were needed to adapt the following MATLAB code to fortran
        !F=Icase==14 | Icase==15; % input TA, (pCO2 or fCO2)
        !if any(F)
        !    PHic(F) = CalculatepHfromTAfCO2(TAc(F)-PengCorrection(F), FCic(F));
        !    TCc(F)  = CalculateTCfromTApH  (TAc(F)-PengCorrection(F), PHic(F));
        !end
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-3'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==14 .or. Icase==15)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            allocate(TAc_minus_Peng(sum(INDEX_ARRAY)), TCc_Interface (sum(INDEX_ARRAY)), &
                     PHic_Interface(sum(INDEX_ARRAY)), FCic_Interface(sum(INDEX_ARRAY)))

            TAc_minus_Peng = TAc (INDEX_ARRAY) - PengCorrection(INDEX_ARRAY)
            FCic_Interface = PHic(INDEX_ARRAY)

            call CalculatepHfromTAfCO2 &
                 (TAc_minus_Peng, FCic_Interface, PHic_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            call CalculateTCfromTApH  &
                 (TAc_minus_Peng, PHic_Interface, TCc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            PHic(INDEX_ARRAY) = PHic_Interface
            TCc (INDEX_ARRAY) = TCc

            deallocate(TAc_minus_Peng, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of input TA, (pCO2 or fCO2)
        ! *******************************************************************************


        ! *******************************************************************************
        ! input TC, pH
        ! *******************************************************************************
        !F=Icase==23;
        ! Changes were needed to adapt the following MATLAB code to fortran
        !if any(F)
        !    TAc(F)  = CalculateTAfromTCpH  (TCc(F), PHic(F)) + PengCorrection(F);
        !    FCic(F) = CalculatefCO2fromTCpH(TCc(F), PHic(F));
        !end
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-4'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==23)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            allocate(TAc_interface (sum(INDEX_ARRAY)), TCc_Interface (sum(INDEX_ARRAY)), &
                     PHic_Interface(sum(INDEX_ARRAY)), FCic_Interface(sum(INDEX_ARRAY)))

            PHic_Interface = PHic(INDEX_ARRAY)
            TCc_Interface  = TCc (INDEX_ARRAY)

            call CalculateTAfromTCpH &
                 (TCc_Interface, PHic_Interface, TAc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            TAc_Interface = TAc_Interface + PengCorrection(INDEX_ARRAY)

            call CalculatefCO2fromTCpH &
                 (TCc_Interface, PHic_Interface, FCic_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            TAc (INDEX_ARRAY) = TAc_Interface
            FCic(INDEX_ARRAY) = FCic_Interface

            deallocate(TAc_interface, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of input TC, pH
        ! *******************************************************************************

        ! *******************************************************************************
        ! input TC, (pCO2 or fCO2)
        ! *******************************************************************************

        ! Changes were needed to adapt the following MATLAB code to fortran
        !F=Icase==24 | Icase==25;  % input TC, (pCO2 or fCO2)
        !if any(F)
        !    PHic(F) = CalculatepHfromTCfCO2(TCc(F), FCic(F));
        !    TAc(F)  = CalculateTAfromTCpH  (TCc(F), PHic(F)) + PengCorrection(F);
        !end
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-5'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==24 .or. Icase==25)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            allocate(TAc_interface (sum(INDEX_ARRAY)), TCc_Interface (sum(INDEX_ARRAY)), &
                     PHic_Interface(sum(INDEX_ARRAY)), FCic_Interface(sum(INDEX_ARRAY)))

            FCic_Interface = FCic(INDEX_ARRAY)
            TCc_Interface  = TCc (INDEX_ARRAY)

            call CalculatepHfromTCfCO2 &
                 (TCc_Interface, FCic_Interface, PHic_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            call CalculateTAfromTCpH  &
                 (TCc_Interface, PHic_Interface, TAc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            TAc_Interface = TAc_Interface + PengCorrection

            TAc (INDEX_ARRAY) = TAc_Interface
            PHic(INDEX_ARRAY) = PHic_Interface

            deallocate(TAc_interface, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of input TC, (pCO2 or fCO2)
        ! *******************************************************************************


        ! *******************************************************************************
        ! input pH, (pCO2 or fCO2)
        ! *******************************************************************************

        ! Changes were needed to adapt the following MATLAB code to fortran
        !F=Icase==34 | Icase==35; % input pH, (pCO2 or fCO2)
        !if any(F)
        !    TCc(F)  = CalculateTCfrompHfCO2(PHic(F), FCic(F));
        !    TAc(F)  = CalculateTAfromTCpH  (TCc(F),  PHic(F)) + PengCorrection(F);
        !end
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #7-6'
        endif

        INDEX_ARRAY_MASK = 0

        where(Icase==34 .or. Icase==35)
            INDEX_ARRAY_MASK = 1
        end where

        if (sum(INDEX_ARRAY_MASK) > 0) then
            call GENERATE_INDEX_ARRAY(INDEX_ARRAY_MASK, INDEX_ARRAY)

            allocate(TAc_interface (sum(INDEX_ARRAY)), TCc_Interface (sum(INDEX_ARRAY)), &
                     PHic_Interface(sum(INDEX_ARRAY)), FCic_Interface(sum(INDEX_ARRAY)))

            PHic_Interface = PHic(INDEX_ARRAY)
            FCic_Interface = FCic(INDEX_ARRAY)

            call CalculateTCfrompHfCO2 &
                 (PHic_Interface, FCic_Interface, TCc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            call CalculateTAfromTCpH  &
                 (TCc_Interface , PHic_Interface, TAc_Interface, &
                  ntps, &
                  pHScale, WhichKs, WhoseKSO4, p1, p2, &
                  TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                  RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                  fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                  FugFac, VPFac)

            TAc_Interface = TAc_Interface + PengCorrection

            TAc(INDEX_ARRAY) = TAc_Interface
            TCc(INDEX_ARRAY) = TCc_Interface

            deallocate(TAc_interface, TCc_Interface, PHic_Interface, FCic_Interface)
        endif
        ! *******************************************************************************
        ! End of input pH, (pCO2 or fCO2)
        ! *******************************************************************************

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CO2SYS : Checkpoint #8'
        endif


        ! By now, an fCO2 value is available for each sample.
        ! Generate the associated pCO2 value:
        ! MATLAB CODE : PCic = FCic./FugFac;
        
        if(.not. allocated(PCic)) allocate(PCic(ntps))
        !call ASSIGN_DBL_VECTOR_CONTENT(PCic, FCic / FugFac)        
        PCic = FCic / FugFac

        ! CalculateOtherParamsAtInputConditions:
        ! MATLAB CODE:
        ! [HCO3inp CO3inp BAlkinp...
        !     OHinp PAlkinp...
        !     SiAlkinp Hfreeinp...
        !     HSO4inp HFinp]      = CalculateAlkParts(PHic, TCc);
        ! PAlkinp                 = PAlkinp+PengCorrection;
        ! CO2inp                  = TCc - CO3inp - HCO3inp;
        ! F=true(ntps,1);           % i.e., do for all samples:
        ! Revelleinp              = RevelleFactor(TAc-PengCorrection, TCc);
        ! [OmegaCainp OmegaArinp] = CaSolubility(Sal, TempCi, Pdbari, TCc, PHic);
        ! xCO2dryinp              = PCic./VPFac; % ' this assumes pTot = 1 atm
        call CalculateAlkParts &
             (PHic, TCc, HCO3inp, CO3inp, BAlkinp, OHinp, PAlkinp, SiAlkinp, &
              Hfreeinp, HSO4inp, HFinp, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call ASSIGN_DBL_VECTOR_CONTENT(PAlkinp, PAlkinp + PengCorrection)
        call ASSIGN_DBL_VECTOR_CONTENT(CO2inp, TCc - CO3inp - HCO3inp)

        call ASSIGN_DBL_VECTOR_CONTENT(TAc_minus_Peng, TAc-PengCorrection)

        call RevelleFactor &
             (DEBUG, TAc_minus_Peng, TCc, Revelleinp, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call CaSolubility &
             (Sal, TempCi, Pdbari, TCc, PHic, OmegaCainp, OmegaArinp, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call ASSIGN_DBL_VECTOR_CONTENT(xCO2dryinp, PCic / VPFac) ! ' this assumes pTot = 1 atm

        ! Just for reference, convert pH at input conditions to the other scales, too.
        ! MATLAB CODE : [pHicT pHicS pHicF pHicN]=FindpHOnAllScales(PHic);
        call FindpHOnAllScales &
             (PHic, pHicT, pHicS, pHicF, pHicN, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        ! Merge the Ks at input into an array. Ks at output will be glued to this later.
        ! MATLAB CODE: KIVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi];
        allocate(KIVEC(size(K0), 13))

        KIVEC(:,1)  = K0(:)
        KIVEC(:,2)  = K1(:)
        KIVEC(:,3)  = K2(:)
        KIVEC(:,4)  = (-log10(K1))
        KIVEC(:,5)  = (-log10(K2))
        KIVEC(:,6)  = KW
        KIVEC(:,7)  = KB
        KIVEC(:,8)  = KF
        KIVEC(:,9)  = KS
        KIVEC(:,10) = KP1
        KIVEC(:,11) = KP2
        KIVEC(:,12) = KP3
        KIVEC(:,13) = KSi

        ! Calculate the constants for all samples at output conditions
        ! MATLAB CODE: Constants(TempCo,Pdbaro);
        call Constants &
             (TempCo, Pdbaro, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        ! Calculate, for output conditions, using conservative TA and TC, pH, fCO2 and pCO2
        ! MATLAB CODE:
        ! F=true(ntps,1); % i.e., do for all samples:
        ! [PHoc FCoc] = CalculatepHfCO2fromTATC(TAc-PengCorrection, TCc);
        ! PCoc = FCoc./FugFac

        call ASSIGN_DBL_VECTOR_CONTENT(TAc_minus_Peng, TAc-PengCorrection)

        call CalculatepHfCO2fromTATC &
             (DEBUG, TAc_minus_Peng, TCc, PHoc, FCoc, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call ASSIGN_DBL_VECTOR_CONTENT(PCoc, FCoc / FugFac)

        ! Calculate Other Stuff At Output Conditions:
        ! MATLAB CODE:
        ! [HCO3out CO3out BAlkout...
        !    OHout PAlkout...
        !    SiAlkout Hfreeout...
        !    HSO4out HFout]      = CalculateAlkParts(PHoc, TCc);
        !PAlkout                 = PAlkout+PengCorrection;
        !CO2out                  = TCc - CO3out - HCO3out;
        !Revelleout              = RevelleFactor(TAc, TCc);
        ![OmegaCaout OmegaArout] = CaSolubility(Sal, TempCo, Pdbaro, TCc, PHoc);
        !xCO2dryout              = PCoc./VPFac; % ' this assumes pTot = 1 atm

        call CalculateAlkParts &
             (PHoc, TCc, HCO3out, CO3out, BAlkout, OHout, PAlkout, SiAlkout, &
              Hfreeout, HSO4out, HFout,&
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call ASSIGN_DBL_VECTOR_CONTENT(PAlkout, PAlkout + PengCorrection)
        call ASSIGN_DBL_VECTOR_CONTENT(CO2out, TCc - CO3out - HCO3out)

        call ASSIGN_DBL_VECTOR_CONTENT(TAc_minus_Peng, TAc-PengCorrection)

        call RevelleFactor &
             (DEBUG, TAc, TCc, Revelleout, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call CaSolubility &
             (Sal, TempCo, Pdbaro, TCc, PHoc, OmegaCaout, OmegaArout, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call ASSIGN_DBL_VECTOR_CONTENT(xCO2dryout, PCoc / VPFac) ! ' this assumes pTot = 1 atm

        ! Just for reference, convert pH at output conditions to the other scales, too.
        ! MATLAB CODE: [pHocT pHocS pHocF pHocN]=FindpHOnAllScales(PHoc);
        call FindpHOnAllScales &
             (PHoc, pHocT, pHocS, pHocF, pHocN, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)


        !MATLAB CODE: KOVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi];
        allocate(KOVEC(size(K0), 13))

        KOVEC(:,1)  = K0(:)
        KOVEC(:,2)  = K1(:)
        KOVEC(:,3)  = K2(:)
        KOVEC(:,4)  = (-log10(K1))
        KOVEC(:,5)  = (-log10(K2))
        KOVEC(:,6)  = KW
        KOVEC(:,7)  = KB
        KOVEC(:,8)  = KF
        KOVEC(:,9)  = KS
        KOVEC(:,10) = KP1
        KOVEC(:,11) = KP2
        KOVEC(:,12) = KP3
        KOVEC(:,13) = KSi

        !MATLAB CODE: TVEC =[TB TF TS];
        allocate(TVEC(size(K0), 3))

        TVEC(:,1) = TB(:)
        TVEC(:,2) = TF(:)
        TVEC(:,3) = TS(:)

        ! Saving data in array, 81 columns, as many rows as samples input
        ! MATLAB CODE:
        !DATA=[TAc*1e6        TCc*1e6         PHic           PCic*1e6        FCic*1e6...
        !      HCO3inp*1e6    CO3inp*1e6      CO2inp*1e6     BAlkinp*1e6     OHinp*1e6...
        !      PAlkinp*1e6    SiAlkinp*1e6    Hfreeinp*1e6   Revelleinp      OmegaCainp... %%% Multiplied Hfreeinp *1e6, svh20100827
        !      OmegaArinp     xCO2dryinp*1e6  PHoc           PCoc*1e6        FCoc*1e6...
        !      HCO3out*1e6    CO3out*1e6      CO2out*1e6     BAlkout*1e6     OHout*1e6...
        !      PAlkout*1e6    SiAlkout*1e6    Hfreeout*1e6   Revelleout      OmegaCaout... %%% Multiplied Hfreeout *1e6, svh20100827
        !      OmegaArout     xCO2dryout*1e6  pHicT          pHicS           pHicF...
        !      pHicN          pHocT           pHocS          pHocF           pHocN...
        !      TEMPIN         TEMPOUT         PRESIN         PRESOUT         PAR1TYPE...
        !      PAR2TYPE       K1K2CONSTANTS   KSO4CONSTANTS  pHSCALEIN       SAL...
        !      PO4            SI              KIVEC          KOVEC           TVEC*1e6];

        if (allocated(CO2SYS_OUT_DATA)) deallocate(CO2SYS_OUT_DATA)
        allocate(CO2SYS_OUT_DATA(size(TAc), 81))

        CO2SYS_OUT_DATA(:,  1  ) = TAc          (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  2  ) = TCc          (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  3  ) = PHic         (:)
        CO2SYS_OUT_DATA(:,  4  ) = PCic         (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  5  ) = FCic         (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  6  ) = HCO3inp      (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  7  ) = CO3inp       (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  8  ) = CO2inp       (:) * 1.0D6
        CO2SYS_OUT_DATA(:,  9  ) = BAlkinp      (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 10  ) = OHinp        (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 11  ) = PAlkinp      (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 12  ) = SiAlkinp     (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 13  ) = Hfreeinp     (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 14  ) = Revelleinp   (:)
        CO2SYS_OUT_DATA(:, 15  ) = OmegaCainp   (:) ! Multiplied Hfreeinp *1.0D6, svh20100827
        CO2SYS_OUT_DATA(:, 16  ) = OmegaArinp   (:)
        CO2SYS_OUT_DATA(:, 17  ) = xCO2dryinp   (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 18  ) = PHoc         (:)
        CO2SYS_OUT_DATA(:, 19  ) = PCoc         (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 20  ) = FCoc         (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 21  ) = HCO3out      (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 22  ) = CO3out       (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 23  ) = CO2out       (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 24  ) = BAlkout      (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 25  ) = OHout        (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 26  ) = PAlkout      (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 27  ) = SiAlkout     (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 28  ) = Hfreeout     (:) * 1.0D6
        CO2SYS_OUT_DATA(:, 29  ) = Revelleout   (:)
        CO2SYS_OUT_DATA(:, 30  ) = OmegaCaout   (:)  !%% Multiplied Hfreeout *1.0D6, svh20100827
        CO2SYS_OUT_DATA(:, 31  ) = OmegaArout   (:)
        CO2SYS_OUT_DATA(:, 32  ) = xCO2dryout   (:) *1.0D6
        CO2SYS_OUT_DATA(:, 33  ) = pHicT        (:)
        CO2SYS_OUT_DATA(:, 34  ) = pHicS        (:)
        CO2SYS_OUT_DATA(:, 35  ) = pHicF        (:)
        CO2SYS_OUT_DATA(:, 36  ) = pHicN        (:)
        CO2SYS_OUT_DATA(:, 37  ) = pHocT        (:)
        CO2SYS_OUT_DATA(:, 38  ) = pHocS        (:)
        CO2SYS_OUT_DATA(:, 39  ) = pHocF        (:)
        CO2SYS_OUT_DATA(:, 40  ) = pHocN        (:)
        CO2SYS_OUT_DATA(:, 41  ) = TEMPIN       (:)
        CO2SYS_OUT_DATA(:, 42  ) = TEMPOUT      (:)
        CO2SYS_OUT_DATA(:, 43  ) = PRESIN       (:)
        CO2SYS_OUT_DATA(:, 44  ) = PRESOUT      (:)
        CO2SYS_OUT_DATA(:, 45  ) = PAR1TYPE     (:)
        CO2SYS_OUT_DATA(:, 46  ) = PAR2TYPE     (:)
        CO2SYS_OUT_DATA(:, 47  ) = K1K2CONSTANTS(:)
        CO2SYS_OUT_DATA(:, 48  ) = KSO4CONSTANTS(:)
        CO2SYS_OUT_DATA(:, 49  ) = pHSCALEIN    (:)
        CO2SYS_OUT_DATA(:, 50  ) = SAL          (:)
        CO2SYS_OUT_DATA(:, 51  ) = PO4          (:)
        CO2SYS_OUT_DATA(:, 52  ) = SI           (:)
        CO2SYS_OUT_DATA(:, 53:65) = KIVEC       (:,:)
        CO2SYS_OUT_DATA(:, 66:78) = KOVEC       (:,:)
        CO2SYS_OUT_DATA(:, 79:81) = TVEC        (:,:) * 1.0D6

        if (allocated(NICEHEADERS)) deallocate(NICEHEADERS)
        allocate(NICEHEADERS(81))

        NICEHEADERS( 1) = '        01 - TAlk (micro mol/kgSW)'
        NICEHEADERS( 2) = '        02 - TCO2 (micro mol/kgSW)'
        NICEHEADERS( 3) = '                      03 - pHin ()'
        NICEHEADERS( 4) = '           04 - pCO2in (micro atm)'
        NICEHEADERS( 5) = '           05 - fCO2in (micro atm)'
        NICEHEADERS( 6) = '      06 - HCO3in (micro mol/kgSW)'
        NICEHEADERS( 7) = '       07 - CO3in (micro mol/kgSW)'
        NICEHEADERS( 8) = '       08 - CO2in (micro mol/kgSW)'
        NICEHEADERS( 9) = '      09 - BAlkin (micro mol/kgSW)'
        NICEHEADERS(10) = '        10 - OHin (micro mol/kgSW)'
        NICEHEADERS(11) = '      11 - PAlkin (micro mol/kgSW)'
        NICEHEADERS(12) = '     12 - SiAlkin (micro mol/kgSW)'
        NICEHEADERS(13) = '     13 - Hfreein (micro mol/kgSW)'
        NICEHEADERS(14) = '          14 - RevelleFactorin  ()'
        NICEHEADERS(15) = '                 15 - OmegaCain ()'
        NICEHEADERS(16) = '                 16 - OmegaArin ()'
        NICEHEADERS(17) = '                 17 - xCO2in (ppm)'
        NICEHEADERS(18) = '                     18 - pHout ()'
        NICEHEADERS(19) = '          19 - pCO2out (micro atm)'
        NICEHEADERS(20) = '          20 - fCO2out (micro atm)'
        NICEHEADERS(21) = '     21 - HCO3out (micro mol/kgSW)'
        NICEHEADERS(22) = '      22 - CO3out (micro mol/kgSW)'
        NICEHEADERS(23) = '      23 - CO2out (micro mol/kgSW)'
        NICEHEADERS(24) = '     24 - BAlkout (micro mol/kgSW)'
        NICEHEADERS(25) = '       25 - OHout (micro mol/kgSW)'
        NICEHEADERS(26) = '     26 - PAlkout (micro mol/kgSW)'
        NICEHEADERS(27) = '    27 - SiAlkout (micro mol/kgSW)'
        NICEHEADERS(28) = '    28 - Hfreeout (micro mol/kgSW)'
        NICEHEADERS(29) = '          29 - RevelleFactorout ()'
        NICEHEADERS(30) = '                30 - OmegaCaout ()'
        NICEHEADERS(31) = '                31 - OmegaArout ()'
        NICEHEADERS(32) = '                32 - xCO2out (ppm)'
        NICEHEADERS(33) = '              33 - pHin (Total) ()'
        NICEHEADERS(34) = '                34 - pHin (SWS) ()'
        NICEHEADERS(35) = '               35 - pHin (Free) ()'
        NICEHEADERS(36) = '               36 - pHin (NBS ) ()'
        NICEHEADERS(37) = '              37 - pHout(Total) ()'
        NICEHEADERS(38) = '                38 - pHout(SWS) ()'
        NICEHEADERS(39) = '               39 - pHout(Free) ()'
        NICEHEADERS(40) = '               40 - pHout(NBS ) ()'
        NICEHEADERS(41) = '               41 - TEMPIN (Deg C)'
        NICEHEADERS(42) = '              42 - TEMPOUT (Deg C)'
        NICEHEADERS(43) = '                43 - PRESIN (dbar)'
        NICEHEADERS(44) = '               44 - PRESOUT (dbar)'
        NICEHEADERS(45) = '                  45 - PAR1TYPE ()'
        NICEHEADERS(46) = '                  46 - PAR2TYPE ()'
        NICEHEADERS(47) = '             47 - K1K2CONSTANTS ()'
        NICEHEADERS(48) = '             48 - KSO4CONSTANTS ()'
        NICEHEADERS(49) = '                 49 - pHSCALEIN ()'
        NICEHEADERS(50) = '         50 - SAL (micro mol/kgSW)'
        NICEHEADERS(51) = '         51 - PO4 (micro mol/kgSW)'
        NICEHEADERS(52) = '          52 - SI (micro mol/kgSW)'
        NICEHEADERS(53) = '                   53 - K0input ()'
        NICEHEADERS(54) = '                   54 - K1input ()'
        NICEHEADERS(55) = '                   55 - K2input ()'
        NICEHEADERS(56) = '                  56 - pK1input ()'
        NICEHEADERS(57) = '                  57 - pK2input ()'
        NICEHEADERS(58) = '                   58 - KWinput ()'
        NICEHEADERS(59) = '                   59 - KBinput ()'
        NICEHEADERS(60) = '                   60 - KFinput ()'
        NICEHEADERS(61) = '                   61 - KSinput ()'
        NICEHEADERS(62) = '                  62 - KP1input ()'
        NICEHEADERS(63) = '                  63 - KP2input ()'
        NICEHEADERS(64) = '                  64 - KP3input ()'
        NICEHEADERS(65) = '                  65 - KSiinput ()'
        NICEHEADERS(66) = '                  66 - K0output ()'
        NICEHEADERS(67) = '                  67 - K1output ()'
        NICEHEADERS(68) = '                  68 - K2output ()'
        NICEHEADERS(69) = '                 69 - pK1output ()'
        NICEHEADERS(70) = '                 70 - pK2output ()'
        NICEHEADERS(71) = '                  71 - KWoutput ()'
        NICEHEADERS(72) = '                  72 - KBoutput ()'
        NICEHEADERS(73) = '                  73 - KFoutput ()'
        NICEHEADERS(74) = '                  74 - KSoutput ()'
        NICEHEADERS(75) = '                 75 - KP1output ()'
        NICEHEADERS(76) = '                 76 - KP2output ()'
        NICEHEADERS(77) = '                 77 - KP3output ()'
        NICEHEADERS(78) = '                 78 - KSioutput ()'
        NICEHEADERS(79) = '          79 - TB (micro mol/kgSW)'
        NICEHEADERS(80) = '          80 - TF (micro mol/kgSW)'
        NICEHEADERS(81) = '          81 - TS (micro mol/kgSW)'

        !Deallocate all internal arrays
        deallocate(TA              )     !if(allocated(TA              ))
        deallocate(TC              )     !if(allocated(TC              ))
        deallocate(PH              )     !if(allocated(PH              ))
        deallocate(PC              )     !if(allocated(PC              ))
        deallocate(FC              )     !if(allocated(FC              ))
        deallocate(PengCorrection  )     !if(allocated(PengCorrection  ))
        deallocate(TAc             )     !if(allocated(TAc             ))        
        deallocate(TAc_minus_Peng  )     !if(allocated(TAc_minus_Peng  ))
        deallocate(TCc             )     !if(allocated(TCc             ))        
        deallocate(PHic            )     !if(allocated(PHic            ))        
        deallocate(PCic            )     !if(allocated(PCic            ))        
        deallocate(FCic            )     !if(allocated(FCic            ))        
        deallocate(Icase           )     !if(allocated(Icase           ))
        deallocate(INDEX_ARRAY     )     !if(allocated(INDEX_ARRAY     ))
        deallocate(INDEX_ARRAY_MASK)     !if(allocated(INDEX_ARRAY_MASK))
        deallocate(CO2inp          )     !if(allocated(CO2inp          ))
        deallocate(HCO3inp         )     !if(allocated(HCO3inp         ))
        deallocate(CO3inp          )     !if(allocated(CO3inp          ))
        deallocate(BAlkinp         )     !if(allocated(BAlkinp         ))
        deallocate(OHinp           )     !if(allocated(OHinp           ))
        deallocate(PAlkinp         )     !if(allocated(PAlkinp         ))
        deallocate(SiAlkinp        )     !if(allocated(SiAlkinp        ))
        deallocate(Hfreeinp        )     !if(allocated(Hfreeinp        ))
        deallocate(HSO4inp         )     !if(allocated(HSO4inp         ))
        deallocate(HFinp           )     !if(allocated(HFinp           ))
        deallocate(Revelleinp      )     !if(allocated(Revelleinp      ))
        deallocate(OmegaCainp      )     !if(allocated(OmegaCainp      ))
        deallocate(OmegaArinp      )     !if(allocated(OmegaArinp      ))
        deallocate(xCO2dryinp      )     !if(allocated(xCO2dryinp      ))
        deallocate(pHicT           )     !if(allocated(pHicT           ))
        deallocate(pHicS           )     !if(allocated(pHicS           ))
        deallocate(pHicF           )     !if(allocated(pHicF           ))
        deallocate(pHicN           )     !if(allocated(pHicN           ))
        deallocate(KIVEC           )     !if(allocated(KIVEC           ))
        deallocate(KOVEC           )     !if(allocated(KOVEC           ))
        deallocate(TVEC            )     !if(allocated(TVEC            ))
        deallocate(pHoc            )     !if(allocated(pHoc            ))
        deallocate(pHocT           )     !if(allocated(pHocT           ))
        deallocate(pHocS           )     !if(allocated(pHocS           ))
        deallocate(pHocF           )     !if(allocated(pHocF           ))
        deallocate(pHocN           )     !if(allocated(pHocN           ))
        deallocate(FCoc            )     !if(allocated(FCoc            ))
        deallocate(PCoc            )     !if(allocated(PCoc            ))
        deallocate(CO2out          )     !if(allocated(CO2out          ))
        deallocate(HCO3out         )     !if(allocated(HCO3out         ))
        deallocate(CO3out          )     !if(allocated(CO3out          ))
        deallocate(BAlkout         )     !if(allocated(BAlkout         ))
        deallocate(OHout           )     !if(allocated(OHout           ))
        deallocate(PAlkout         )     !if(allocated(PAlkout         ))
        deallocate(SiAlkout        )     !if(allocated(SiAlkout        ))
        deallocate(Hfreeout        )     !if(allocated(Hfreeout        ))
        deallocate(HSO4out         )     !if(allocated(HSO4out         ))
        deallocate(HFout           )     !if(allocated(HFout           ))
        deallocate(Revelleout      )     !if(allocated(Revelleout      ))
        deallocate(OmegaCaout      )     !if(allocated(OmegaCaout      ))
        deallocate(OmegaArout      )     !if(allocated(OmegaArout      ))
        deallocate(xCO2dryout      )     !if(allocated(xCO2dryout      ))


        !Deallocate all former global variables of CO2SYS
         deallocate(pHScale  )      !if(allocated(pHScale  ))
         deallocate(WhichKs  )      !if(allocated(WhichKs  ))
         deallocate(WhoseKSO4)      !if(allocated(WhoseKSO4))
         deallocate(p1       )      !if(allocated(p1       ))
         deallocate(p2       )      !if(allocated(p2       ))
         deallocate(TempCi   )      !if(allocated(TempCi   ))
         deallocate(TempCo   )      !if(allocated(TempCo   ))
         deallocate(Pdbari   )      !if(allocated(Pdbari   ))
         deallocate(Pdbaro   )      !if(allocated(Pdbaro   ))
         deallocate(Sal      )      !if(allocated(Sal      ))
         deallocate(sqrSal   )      !if(allocated(sqrSal   ))
         deallocate(TP       )      !if(allocated(TP       ))
         deallocate(TSi      )      !if(allocated(TSi      ))
         deallocate(TB       )      !if(allocated(TB       ))
         deallocate(TF       )      !if(allocated(TF       ))
         deallocate(TS       )      !if(allocated(TS       ))
         deallocate(TempK    )      !if(allocated(TempK    ))
         deallocate(RT       )      !if(allocated(RT       ))
         deallocate(logTempK )      !if(allocated(logTempK ))
         deallocate(Pbar     )      !if(allocated(Pbar     ))
         deallocate(TempK100 )      !if(allocated(TempK100 ))
         deallocate(lnK0     )      !if(allocated(lnK0     ))
         deallocate(K0       )      !if(allocated(K0       ))
         deallocate(KS       )      !if(allocated(KS       ))
         deallocate(fH       )      !if(allocated(fH       ))
         deallocate(K1       )      !if(allocated(K1       ))
         deallocate(K2       )      !if(allocated(K2       ))
         deallocate(KW       )      !if(allocated(KW       ))
         deallocate(KB       )      !if(allocated(KB       ))
         deallocate(KF       )      !if(allocated(KF       ))
         deallocate(KP1      )      !if(allocated(KP1      ))
         deallocate(KP2      )      !if(allocated(KP2      ))
         deallocate(KP3      )      !if(allocated(KP3      ))
         deallocate(KSi      )      !if(allocated(KSi      ))
         deallocate(FugFac   )      !if(allocated(FugFac   ))
         deallocate(VPFac    )      !if(allocated(VPFac    ))
        !End of deallocate all former global variables of CO2SYS
        
    end subroutine CO2SYS

    ! **************************************************************************************
    ! Subroutine that calculates the constants
    ! **************************************************************************************

    ! Historical notes from original CO2SYS implemented in QBASIC
    ! SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
    ! Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
    ! Outputs: K0, K(), T(), fH, FugFac, VPFac
    ! This finds the Constants of the CO2 system in seawater or freshwater,
    ! corrects them for pressure, and reports them on the chosen pH scale.
    ! The process is as follows: the Constants (except KS, KF which stay on the
    ! free scale - these are only corrected for pressure) are
    !        1) evaluated as they are given in the literature
    !        2) converted to the SWS scale in mol/kg-SW or to the NBS scale
    !        3) corrected for pressure
    !        4) converted to the SWS pH scale in mol/kg-SW
    !        5) converted to the chosen pH scale
    !
    !        PROGRAMMER'S NOTE: all logs are log base e
    !        PROGRAMMER'S NOTE: all Constants are converted to the pH scale
    !                pHScale% (the chosen one) in units of mol/kg-SW
    !                except KS and KF are on the free scale
    !                and KW is in units of (mol/kg-SW)^2
    subroutine Constants &
               (TempC, Pdbar, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! Argument list
        real(kind = DBL_PREC), intent(in), dimension(:) :: TempC
        real(kind = DBL_PREC), intent(in), dimension(:) :: Pdbar
        ! End of argument list

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS


        ! Auxillary variables
        real(kind = DBL_PREC), allocatable, dimension(:) :: IonS
        real(kind = DBL_PREC), allocatable, dimension(:) :: pKS
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKS

        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKF
        real(kind = DBL_PREC), allocatable, dimension(:) :: SWStoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT

        real(kind = DBL_PREC), allocatable, dimension(:) :: logKB
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKBtop
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKB

        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKW

        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP1
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP2
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP3
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKSi

        real(kind = DBL_PREC), allocatable, dimension(:) :: logK1
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnK1
        real(kind = DBL_PREC), allocatable, dimension(:) :: pK1
        real(kind = DBL_PREC), allocatable, dimension(:) :: logK2
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnK2
        real(kind = DBL_PREC), allocatable, dimension(:) :: pK2

        real(kind = DBL_PREC), allocatable, dimension(:) :: F1
        real(kind = DBL_PREC), allocatable, dimension(:) :: F2
        real(kind = DBL_PREC), allocatable, dimension(:) :: PK1_0
        real(kind = DBL_PREC), allocatable, dimension(:) :: A_1
        real(kind = DBL_PREC), allocatable, dimension(:) :: B_1
        real(kind = DBL_PREC), allocatable, dimension(:) :: C_1
        real(kind = DBL_PREC), allocatable, dimension(:) :: PK2_0
        real(kind = DBL_PREC), allocatable, dimension(:) :: A_2
        real(kind = DBL_PREC), allocatable, dimension(:) :: B_2
        real(kind = DBL_PREC), allocatable, dimension(:) :: C_2

        real(kind = DBL_PREC), allocatable, dimension(:) :: PK10
        real(kind = DBL_PREC), allocatable, dimension(:) :: A1
        real(kind = DBL_PREC), allocatable, dimension(:) :: B1
        real(kind = DBL_PREC), allocatable, dimension(:) :: C1
        real(kind = DBL_PREC), allocatable, dimension(:) :: PK20
        real(kind = DBL_PREC), allocatable, dimension(:) :: A2
        real(kind = DBL_PREC), allocatable, dimension(:) :: B2
        real(kind = DBL_PREC), allocatable, dimension(:) :: C2

        real(kind = DBL_PREC), allocatable, dimension(:) :: deltaV
        real(kind = DBL_PREC), allocatable, dimension(:) :: Kappa
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnK1fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnK2fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKBfac
        
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKWfac

        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKFfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKSfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP1fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP2fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKP3fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKSifac
        real(kind = DBL_PREC), allocatable, dimension(:) :: K1fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KWfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KBfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KFfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSfac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3fac
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSifac
        
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHfactor
        real(kind = DBL_PREC), allocatable, dimension(:) :: Delta
        real(kind = DBL_PREC), allocatable, dimension(:) :: b
        real(kind = DBL_PREC) :: P1atm
        real(kind = DBL_PREC), allocatable, dimension(:) :: VPWP
        real(kind = DBL_PREC), allocatable, dimension(:) :: VPCorrWP
        real(kind = DBL_PREC), allocatable, dimension(:) :: VPSWWP
        ! End of auxillary variables

        TempK   = TempC + 273.15
        RT      = RGasConstant * TempK
        
        !call ASSIGN_DBL_VECTOR_CONTENT(logTempK, log(TempK))
        !if(.not. allocated(logTempK)) allocate(logTempK(ntps))
        logTempK = log(TempK)
        
        Pbar    = Pdbar / 10
        
        ! Generate empty vectors for holding results
        TB = 0.0D0
        TF = 0.0D0
        TS = 0.0D0

        if (.not.(allocated(lnKS))) then
            allocate(lnKS(ntps))
            lnKS = 0.0D0
        endif


                                  
        ! Calculate TB - Total Borate:

        ! Pure water case
        where(WhichKs==8)
            TB = 0.0D0
        end where

        where(WhichKs==6 .or. WhichKs==7)
            TB = 0.00041060D0 * (Sal / 35.0D0) ! in mol/kg-SW
            ! this is .00001173.*Sali
            ! this is about 1% lower than Uppstrom's value
            ! Culkin, F., in Chemical Oceanography, ed. Riley and Skirrow, 1965:
            ! GEOSECS references this, but this value is not explicitly given here
        end where

        ! All other cases
        where ((WhichKs/=6 .and. WhichKs /= 7 .and. WhichKs/=8) .and. &
               (WhoseKSO4==1 .or. WhoseKSO4==2))    !If user opted for Uppstrom's values
	        ! Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	        ! this is .000416.*Sali./35. = .0000119.*Sali
	        TB =  0.0004157D0 * (Sal/35.0D0) ! in mol/kg-SW
        end where

        where ((WhichKs/=6 .and. WhichKs /= 7 .and. WhichKs/=8) .and. &
               (WhoseKSO4==3 .or. WhoseKSO4==4))    !If user opted for the new Lee values
		    ! Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
	 	    ! Geochimica Et Cosmochimica Acta 74 (6): 1801–1811.
		    TB =  0.0004326D0 * (Sal/35.0D0) ! in mol/kg-SW
        end where
        ! End of calculate TB - Total Borate

        ! CalculateTF;
        ! Riley, J. P., Deep-Sea Research 12:219-220, 1965:
        ! this is .000068.*Sali./35. = .00000195.*Sali
        TF = (0.000067D0 / 18.998D0) * (Sal /1.80655D0) ! in mol/kg-SW

        ! CalculateTS ;
        ! Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
        ! this is .02824.*Sali./35. = .0008067.*Sali
        TS = (0.14D0 / 96.062D0) * (Sal / 1.80655D0)    ! in mol/kg-SW

        ! CalculateK0
        ! Weiss, R. F., Marine Chemistry 2:203-215, 1974.
        TempK100 = TempK / 100.0D0

        lnK0 = -60.2409D0 + 93.4517D0 / TempK100 + 23.3585D0 * log(TempK100) + &
               Sal * (0.023517D0 - 0.023656D0 * TempK100 + 0.0047036D0 * TempK100 * TempK100)

        K0 = exp(lnK0)                  ! this is in mol/kg-SW/atm

        ! CalculateIonS:
        ! This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
        if (.not. allocated(IonS)) allocate(IonS (ntps))
        IonS = 19.924D0 * Sal / (1000.0D0 - (1.005D0 * Sal))

        ! CalculateKS:
        if (.not. allocated(pKS)) allocate(pKS (ntps))
        !if (.not. allocated(KS)) allocate(KS  (ntps))

        pKS = 0.0D0
        KS  = 0.0D0
        


        where(WhoseKSO4==1 .or. WhoseKSO4==3)
            ! Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
            ! The goodness of fit is .021.
            ! It was given in mol/kg-H2O. I convert it to mol/kg-SW.
            ! TYPO on p. 121: the constant e9 should be e8.
            ! This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
            lnKS = -4276.1D0 / TempK + 141.328D0 - 23.093D0 * logTempK + &
                   (-13856.0D0 / TempK + 324.57D0  - 47.986D0 *logTempK) * sqrt(IonS) + &
                      ( 35474.0D0 / TempK - 771.54D0 + 114.723D0 *logTempK) * IonS + &
                       (-2698.0D0 / TempK) * sqrt(IonS) * IonS + (1776 / TempK) * IonS * IonS

	        KS = exp(lnKS) * (1.0D0 - (0.001005D0 * Sal))
        end where

        where(WhoseKSO4==1 .or. WhoseKSO4==3)
            ! Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
            ! KS was found by titrations with a hydrogen electrode
            ! of artificial seawater containing sulfate (but without F)
            ! at 3 salinities from 20 to 45 and artificial seawater NOT
            ! containing sulfate (nor F) at 16 salinities from 15 to 45,
            ! both at temperatures from 5 to 40 deg C.
            ! KS is on the Free pH scale (inherently so).
            ! It was given in mol/kg-H2O. I convert it to mol/kg-SW.
            ! He finds log(beta) which = my pKS;
            ! his beta is an association constant.
            ! The rms error is .0021 in pKS, or about .5% in KS.
            ! This is equation 20 on p. 33:
            pKS = (647.59D0 / TempK) - 6.3451D0 + (0.019085D0 * TempK) - (0.5208D0 * sqrt(IonS))

            ! this is on the free pH scale in mol/kg-H2O convert to mol/kg-SW
            KS = 10.0D0 ** ((-pKS) * (1 - (0.001005D0 * Sal)))
        end where

        ! CalculateKF:
        if(.not.(allocated(lnKF)))      allocate(lnKF(ntps))
        if(.not.(allocated(SWStoTOT)))  allocate(SWStoTOT(ntps))
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))

        ! Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
        lnKF = (1590.2D0 / TempK) - 12.641D0 + (1.525D0 * (IonS**0.5D0))

        ! this is on the free pH scale in mol/kg-H2O convert to mol/kg-SW
        KF = exp(lnKF) * (1.0D0 - (0.001005D0 * Sal))

        ! Another expression exists for KF: Perez and Fraga 1987.
        ! Not used here since ill defined for low salinity.
        ! (to be used for S: 10-40, T: 9-33)
        ! Nonetheless, P&F87 might actually be better than the fit of D&R79 above,
        ! which is based on only three salinities: [0 26.7 34.6]
        ! lnKF = 874./TempK - 9.68 + 0.111.*Sal.^0.5;
        ! KF   = exp(lnKF);                   % this is on the free pH scale in mol/kg-SW

        ! CalculatepHScaleConversionFactors:
        !       These are NOT pressure-corrected.
        SWStoTOT  = (1.0D0 + (TS / KS)) / (1.0D0 + (TS / KS) + (TF / KF))
        FREEtoTOT =  1.0D0 + (TS / KS)

        ! CalculatefH

        ! Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
        where(WhichKs==8)
            ! this shouldn't occur in the program for this case
            fH = 1
        end where

        where(WhichKs==7)
            fH = 1.29D0 - (0.00204D0 * TempK) + ((0.00046D0 - 0.00000148D0 * TempK) * Sal * Sal)
            ! Peng et al, Tellus 39B:439-458, 1987:
            ! They reference the GEOSECS report, but round the value
            ! given there off so that it is about .008 (1%) lower. It
            ! doesn't agree with the check value they give on p. 456.
        end where

        where((WhichKs.ne.7) .and. (WhichKs.ne.8))
            fH = 1.2948D0 - (0.002036D0 * TempK) + &
                 ((0.0004607D0 - 0.000001475D0 * TempK) * Sal * Sal)
            ! Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition v. 3, 1982 (p. 80)
        end where

        ! CalculateKB:
        if(.not.(allocated(logKB)))   allocate(logKB(ntps))
        if(.not.(allocated(lnKBtop))) allocate(lnKBtop(ntps))
        if(.not.(allocated(lnKB)))    allocate(lnKB(ntps))

        ! Pure water case
        where(WhichKs==8)
            KB = 0.0D0
        end where

        where((WhichKs==6) .or. (WhichKs==7));
            ! This is for GEOSECS and Peng et al. Lyman, John, UCLA Thesis, 1957
            ! fit by Li et al, JGR 74:5507-5525, 1969:
            logKB = -9.26D0 + (0.00886D0 * Sal) + (0.01D0 * TempC)

            ! this is on the NBS scale convert to the SWS scale
            KB = (10**logKB) / fH
        end where

        where((WhichKs.ne.6) .and. (WhichKs.ne.7) .and. (WhichKs.ne.8))
            ! Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
            lnKBtop = (-8966.9D0 - 2890.53D0 * sqrSal) - (77.942D0 * Sal) + &
                      (1.728D0 * sqrSal * Sal) -( 0.0996D0 * (Sal**2))

            lnKB = (lnKBtop / TempK) + (148.0248D0 + 137.1942D0 * sqrSal) + &
                   (1.62142D0 * Sal) + (-24.4344D0 - 25.085D0 * sqrSal - 0.2474D0 * Sal) * &
                   (logTempK + 0.053105D0 * sqrSal * TempK)

            ! This is on the total pH scale in mol/kg-SW convert to SWS pH scale
            KB = exp(lnKB) / SWStoTOT;
        end where


        ! CalculateKW
        if(.not.(allocated(lnKW))) allocate(lnKW(ntps))

        where(WhichKs==7)
            ! Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            lnKW = 148.9802D0 - (13847.26D0 / TempK) - (23.6521D0 * logTempK) + &
                ((-79.2447D0 + (3298.72D0 / TempK) + (12.0408D0 * logTempK)) * sqrSal) - &
                (0.019813D0 * Sal)
        end where

        where(WhichKs==8)
            ! Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            ! refit data of Harned and Owen, The Physical Chemistry of
            ! Electrolyte Solutions, 1958
            lnKW = 148.9802D0 - (13847.26D0 / TempK) - (23.6521D0 * logTempK)
        end where

        where((WhichKs.ne.6) .and. (WhichKs.ne.7) .and. (WhichKs.ne.8))
            ! Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
            ! his check value of 1.6 umol/kg-SW should be 6.2
            lnKW = 148.9802D0 - (13847.26D0 / TempK) - (23.6521D0 * logTempK) + &
                ((-5.977D0 + (118.67D0 / TempK) + (1.0495D0 * logTempK)) * sqrSal) - &
                (0.01615D0 * Sal)
        end where

        ! This is on the SWS pH scale in (mol/kg-SW)^2
        KW = exp(lnKW)

        where(WhichKs==6)
            ! GEOSECS doesn't include OH effects
            KW = 0.0D0
        end where

        ! Calculate KP1 KP2 KP3 KSi
        if(.not.(allocated(lnKP1))) allocate(lnKP1(ntps))
        if(.not.(allocated(lnKP2))) allocate(lnKP2(ntps))
        if(.not.(allocated(lnKP3))) allocate(lnKP3(ntps))
        if(.not.(allocated(lnKSi))) allocate(lnKSi(ntps))

        where(WhichKs==7)
            KP1 = 0.02D0
            ! Peng et al don't include the contribution from this term,
            ! but it is so small it doesn't contribute. It needs to be
            ! kept so that the routines work ok.
            ! KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
            ! Limnology and Oceanography 12:243-252, 1967:
            ! these are only for sals 33 to 36 and are on the NBS scale

            ! this is on the NBS scale convert to SWS scale
            KP2 = exp(-9.039D0 - (1450D0 / TempK)) / fH

            ! this is on the NBS scale convert to SWS scale
            KP3 = exp(4.466D0 - (7276D0 / TempK)) / fH

            ! Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
            ! The Chemical Society (London), Special Publ. 17:751, 1964:
            ! this is on the NBS scale convert to SWS scale
            KSi = 0.0000000004D0 / fH
        end where

        where((WhichKs==6) .or. (WhichKs==8))
            KP1 = 0.0D0
            KP2 = 0.0D0
            KP3 = 0.0D0
            KSi = 0.0D0

            ! Neither the GEOSECS choice nor the freshwater choice
            ! include contributions from phosphate or silicate.
        end where

        where((WhichKs.ne.6) .and. (WhichKs.ne.7) .and. (WhichKs.ne.8))
            ! Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
            ! KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
            ! KSi was given on the SWS pH scale in molal units.
            lnKP1 = (-4576.752D0 / TempK) + 115.54D0 - (18.453D0 * logTempK) + &
                 (((-106.736D0 / TempK) + 0.69171D0) * sqrSal) + &
                 (((-0.65643D0 / TempK) - 0.01844D0) * Sal)

            KP1 = exp(lnKP1)

            lnKP2 = (-8814.715D0 / TempK) + 172.1033D0 - (27.927D0 * logTempK) + &
                ((-160.34D0 / TempK) + 1.3566D0) * sqrSal + &
                (((0.37335D0 / TempK) - 0.05778D0) * Sal)

            KP2 = exp(lnKP2)

            lnKP3 = (-3070.75D0 /TempK) - 18.126D0 + &
                (((17.27039D0 / TempK) + 2.81197D0) * sqrSal) + &
                (((-44.99486D0 / TempK) - 0.09984D0) * Sal)

            KP3 = exp(lnKP3)

            lnKSi = (-8904.2D0 /TempK) + 117.4D0 - (19.334D0 * logTempK) + &
                    (((-458.79D0 / TempK) + 3.5913D0) * sqrt(IonS)) + &
                    (((188.74D0 / TempK) - 1.5998D0) * IonS) +  &
                    (((-12.1652D0 / TempK) + 0.07871D0) * IonS * IonS)

            ! this is on the SWS pH scale in mol/kg-H2O convert to mol/kg-SW
            KSi = exp(lnKSi) * (1.0D0 - (0.001005D0 * Sal))
        end where

        if(.not.(allocated(logK1))) allocate(logK1(ntps))
        if(.not.(allocated(lnK1)))  allocate(lnK1(ntps))
        if(.not.(allocated(pK1)))   allocate(pK1(ntps))
        if(.not.(allocated(logK2))) allocate(logK2(ntps))
        if(.not.(allocated(lnK2)))  allocate(lnK2(ntps))
        if(.not.(allocated(pK2)))   allocate(pK2(ntps))
        if(.not.(allocated(F1)))    allocate(F1(ntps))
        if(.not.(allocated(F2)))    allocate(F2 (ntps))
        if(.not.(allocated(PK1_0))) allocate(PK1_0(ntps))
        if(.not.(allocated(A_1)))   allocate(A_1(ntps))
        if(.not.(allocated(B_1)))   allocate(B_1(ntps))
        if(.not.(allocated(C_1)))   allocate(C_1(ntps))
        if(.not.(allocated(PK2_0))) allocate(PK2_0(ntps))
        if(.not.(allocated(A_2)))   allocate(A_2(ntps))
        if(.not.(allocated(B_2)))   allocate(B_2(ntps))
        if(.not.(allocated(C_2)))   allocate(C_2(ntps))

!        if(.not.(allocated(PK1_0))) allocate(PK1_0(ntps)) corrected by Petras: PK1_0 -> PK10
        if(.not.(allocated(PK10))) allocate(PK10(ntps))
        if(.not.(allocated(A1)))    allocate(A1(ntps))
        if(.not.(allocated(B1)))    allocate(B1(ntps))
        if(.not.(allocated(C1)))    allocate(C1(ntps))
        if(.not.(allocated(PK20)))  allocate(PK20(ntps))
        if(.not.(allocated(A2)))    allocate(A2(ntps))
        if(.not.(allocated(B2)))    allocate(B2(ntps))
        if(.not.(allocated(C2)))    allocate(C2(ntps))
        

        where(WhichKs==1)
            ! ROY et al, Marine Chemistry, 44:249-267, 1993
            ! (see also: Erratum, Marine Chemistry 45:337, 1994
            ! and Erratum, Marine Chemistry 52:183, 1996)
            ! Typo: in the abstract on p. 249: in the eq. for lnK1* the
            ! last term should have S raised to the power 1.5.
            ! They claim standard deviations (p. 254) of the fits as
            ! .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            ! They also claim (p. 258) 2s precisions of .004 in pK1 and
            ! .006 in pK2. These are consistent, but Andrew Dickson
            ! (personal communication) obtained an rms deviation of about
            ! .004 in pK1 and .003 in pK2. This would be a 2s precision
            ! of about 2% in K1 and 1.5% in K2.
            ! T:  0-45  S:  5-45. Total Scale. Artificial sewater.
            ! This is eq. 29 on p. 254 and what they use in their abstract:

            lnK1 = 2.83655D0 - (2307.1266D0 / TempK) - (1.5529413D0 * logTempK) + &
                ((-0.20760841D0 - (4.0484D0 / TempK)) * sqrSal) + (0.08468345D0 * Sal) - &
                (0.00654208D0 * sqrSal * Sal)

            ! this is on the total pH scale in mol/kg-H2O, convert to mol/kg-SW, %
            K1 = exp(lnK1) * ((1.0D0 - (0.001005D0 * Sal)) / SWStoTOT)

            ! This is eq. 30 on p. 254 and what they use in their abstract:
            lnK2 = -9.226508D0 - (3351.6106 / TempK) - (0.2005743 * logTempK) + &
                ((-0.106901773 - (23.9722 / TempK)) * sqrSal) + (0.1130822 * Sal) - &
                (0.00846934 * sqrSal * Sal)

            ! this is on the total pH scale in mol/kg-H2O, convert to mol/kg-SW, convert to SWS pH scale
            K2 = exp(lnK2) *((1 - (0.001005 * Sal)) / SWStoTOT)
        end where

        where(WhichKs==2)
            ! GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
            ! The 2s precision in pK1 is .011, or 2.5% in K1.
            ! The 2s precision in pK2 is .02, or 4.5% in K2.
            ! This is in Table 5 on p. 1652 and what they use in the abstract:
            pK1 = (812.27D0 / TempK) + 3.356D0 - (0.00171D0 * Sal *logTempK) + &
                (0.000091D0 *Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K1 = 10.0D0 ** (-pK1)

            ! This is in Table 5 on p. 1652 and what they use in the abstract:
            pK2 = (1450.87D0 / TempK) + 4.604D0 - (0.00385D0 * Sal * logTempK) + &
                (0.000182D0 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K2 = 10.0D0 ** (-pK2)
        end where
        
         

        where(WhichKs==3)
            ! HANSSON refit BY DICKSON AND MILLERO
            ! Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            ! (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            ! refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            ! and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            ! on the SWS pH scale in mol/kg-SW.
            ! Hansson gave his results on the Total scale (he called it
            ! the seawater scale) and in mol/kg-SW.
            ! Typo in DM on p. 1739 in Table 4: the equation for pK2*
            ! for Hansson should have a .000132 *S^2
            ! instead of a .000116 *S^2.
            ! The 2s precision in pK1 is .013, or 3% in K1.
            ! The 2s precision in pK2 is .017, or 4.1% in K2.

            ! This is from Table 4 on p. 1739.
            pK1 = (851.4D0 / TempK) + 3.237D0 - (0.0106D0 * Sal) + (0.000105D0 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K1 = 10.0D0 ** (-pK1)

            ! This is from Table 4 on p. 1739.
            pK2 = (-3885.4D0 / TempK) + 125.844D0 - (18.141D0 * logTempK) - &
                (0.0192D0 * Sal) + (0.000132D0 * Sal * Sal)

            ! This is on the SWS pH scale in mol/kg-SW
            K2 = 10.0D0 ** (-pK2)
        end where

        where(WhichKs==4)
            ! MEHRBACH refit BY DICKSON AND MILLERO
            ! Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            ! (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            ! refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            ! on the SWS pH scale in mol/kg-SW.
            ! Mehrbach et al gave results on the NBS scale.
            ! The 2s precision in pK1 is .011, or 2.6% in K1.
            ! The 2s precision in pK2 is .020, or 4.6% in K2.
	        ! Valid for salinity 20-40.

            ! This is in Table 4 on p. 1739.
            pK1 = (3670.7 / TempK) - 62.008 + (9.7944 * logTempK) - &
                (0.0118 * Sal) + (0.000116 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K1 = 10.0D0 ** (-pK1)

            ! This is in Table 4 on p. 1739.
            pK2 = (1394.7 / TempK) + 4.777 - (0.0184 * Sal) + (0.000118 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K2 = 10.0D0 ** (-pK2)
        end where

        where(WhichKs==5)
            ! HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
            ! Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            ! (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            ! refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            ! Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            ! and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            ! on the SWS pH scale in mol/kg-SW.
            ! Typo in DM on p. 1740 in Table 5: the second equation
            ! should be pK2* =, not pK1* =.
            ! The 2s precision in pK1 is .017, or 4% in K1.
            ! The 2s precision in pK2 is .026, or 6% in K2.
	        ! Valid for salinity 20-40.

            ! This is in Table 5 on p. 1740.
            pK1 = (845.0D0 / TempK) + 3.248D0 - (0.0098D0 * Sal) + (0.000087D0 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K1 = 10.0D0**(-pK1)

            ! This is in Table 5 on p. 1740.
            pK2 = (1377.3D0 / TempK) + 4.824D0 - (0.0185D0 * Sal) + (0.000122D0 * Sal * Sal)

            ! this is on the SWS pH scale in mol/kg-SW
            K2 = 10.0D0**(-pK2)
        end where

        where((WhichKs==6) .or. (WhichKs==7))
            ! GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
            ! Limnology and Oceanography, 18(6):897-907, 1973.
	        ! I.e., these are the original Mehrbach dissociation constants.
            ! The 2s precision in pK1 is .005, or 1.2% in K1.
            ! The 2s precision in pK2 is .008, or 2% in K2.
            pK1 = - 13.7201D0 + (0.031334D0 * TempK) + (3235.76D0 / TempK) + &
                (1.3D-5 * Sal * TempK) - (0.1032D0 * (Sal**0.5D0))

            ! this is on the NBS scale, convert to SWS scale
            K1 = (10.0D0 ** (-pK1)) / fH

            ! pK2 is not defined for Sal=0, since log10(0)=-inf
            pK2 = 5371.9645D0 + (1.671221D0 * TempK) + (0.22913D0 * Sal) + &
                (18.3802D0 * log10(Sal)) - (128375.28D0 / TempK) - &
                (2194.3055D0 * log10(TempK)) - (8.0944D-4 * Sal * TempK) - &
                (5617.11D0 * (log10(Sal) / TempK)) + (2.136D0 * (Sal /TempK))

            ! this is on the NBS scale, convert to SWS scale
            K2 = (10.0D0 ** (-pK2)) / fH
        end where
        
         

        where(WhichKs==8)
	        ! PURE WATER CASE
            ! Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            ! K1 from refit data from Harned and Davis,
            ! J American Chemical Society, 65:2030-2037, 1943.
            ! K2 from refit data from Harned and Scholes,
            ! J American Chemical Society, 43:1706-1709, 1941.
	        ! This is only to be used for Sal=0 water (note the absence of S in the below formulations)
            ! These are the thermodynamic Constants:
            lnK1 = 290.9097D0 - (14554.21D0 / TempK) - (45.0575D0 * logTempK)
            K1 = exp(lnK1)

            lnK2 = 207.6548D0 - (11843.79D0 / TempK) - (33.6485D0 * logTempK)
            K2 = exp(lnK2)
        end where

        where(WhichKs==9)
            ! From Cai and Wang 1998, for estuarine use.
	        ! Data used in this work is from:
	        ! K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
	        ! K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	        ! Sigma of residuals between fits and above data: Â±0.015, +0.040 for K1 and K2, respectively.
	        ! Sal 0-40, Temp 0.2-30
            ! Limnol. Oceanogr. 43(4) (1998) 657-668
	        ! On the NBS scale
	        ! Their check values for F1 don't work out, not sure if this was correctly published...
	        F1 = (200.1D0 / TempK) + 0.3220

	        pK1 = (3404.71D0 /TempK) + (0.032786D0 * TempK) - 14.8435D0 - &
                  (0.071692D0 * F1 * (Sal**0.5D0)) + (0.0021487D0 * Sal)

            ! this is on the NBS scale, convert to SWS scale (uncertain at low Sal due to junction potential)
            K1  = (10.0D0**(-pK1)) / fH

	        F2 = (-129.24D0 / TempK) + 1.4381D0

	        pK2 = (2902.39D0 / TempK) + (0.02379D0 * TempK) - 6.4980D0 - &
                  (0.3191D0 * F2 * (Sal**0.5D0)) + (0.0198D0 * Sal)

            ! this is on the NBS scale, convert to SWS scale (uncertain at low Sal due to junction potential)
            K2 = (10.0D0**(-pK2)) / fH
        end where

        where(WhichKs==10)
            ! From Lueker, Dickson, Keeling, 2000
        	! This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work.
            ! Mar. Chem. 70 (2000) 105-119

            ! Total scale and kg-sw
            pK1 = (3633.86D0 / TempK) - 61.2172D0 + (9.6777D0 * log(TempK)) - &
                (0.011555D0 * Sal) + (0.0001152D0 * Sal * Sal)

            ! this is on the total pH scale in mol/kg-SW, convert to SWS pH scale
	        K1  = (10.0D0 ** (-pK1)) / SWStoTOT

            pK2 = (471.78D0 / TempK) + 25.929D0 -(3.16967D0 * log(TempK)) - &
                (0.01781D0 * Sal) + (0.0001122D0 * Sal * Sal)

            ! this is on the total pH scale in mol/kg-SW, convert to SWS pH scale
	        K2  = (10.0D0 ** (-pK2)) /SWStoTOT
        end where

        where(WhichKs==11)
	        ! Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
  	        ! sigma for pK1 is reported to be 0.0056
	        ! sigma for pK2 is reported to be 0.010
	        ! This is from the abstract and pages 2536-2537
            pK1 = -43.6977D0 - (0.0129037D0 * Sal) + (1.364D-4 * Sal * Sal) + &
                (2885.378D0 / TempK) + (7.045159D0 * log(TempK))

            pK2 = -452.0940D0 + (13.142162D0 * Sal) - (8.101D-4 * Sal * Sal) + &
                (21263.61D0 / TempK) + (68.483143D0 * log(TempK)) + &
                (((-581.4428D0 * Sal) + (0.259601D0 * Sal * Sal)) / TempK) - &
                (1.967035D0 * Sal *log(TempK))

            ! this is on the SWS pH scale in mol/kg-SW
	        K1 = 10.0D0 ** (-pK1)

            ! this is on the SWS pH scale in mol/kg-SW
	        K2 = 10.0D0 ** (-pK2)
        end where

        where(WhichKs==12)
   	        ! Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	        ! Calculated from overdetermined WOCE-era field measurements
	        ! sigma for pK1 is reported to be 0.005
	        ! sigma for pK2 is reported to be 0.008
	        ! This is from page 1715
            pK1 =  6.359D0 - (0.00664D0 * Sal) - (0.01322D0 * TempC) + (4.989D-5 * TempC * TempC)
            pK2 =  9.867D0 - (0.01314D0 * Sal) - (0.01904D0 * TempC) + (2.448D-5 * TempC * TempC)

            ! this is on the SWS pH scale in mol/kg-SW
	        K1 = 10.0D0**(-pK1)

            ! this is on the SWS pH scale in mol/kg-SW
	        K2 = 10.0D0**(-pK2)
        end where

        where(WhichKs==13)
            ! From Millero 2006 work on pK1 and pK2 from titrations
	        ! Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
            ! S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	        pK1_0 = -126.34048D0 + (6320.813D0 / TempK) + (19.568224D0 * log(TempK))
	        A_1   = (13.4191D0 * (Sal**0.5D0)) + (0.0331D0 * Sal) - (5.33D-5 * Sal * Sal)
	        B_1   = (-530.123D0 * (Sal ** 0.5D0)) - (6.103D0 * Sal)
	        C_1   = -2.06950D0 * (Sal ** 0.5D0)

            ! pK1 sigma = 0.0054
	        pK1 = A_1 + (B_1 /TempK) + (C_1 * log(TempK)) + pK1_0

	        K1 = 10.0D0**(-pK1)

	        pK2_0 = -90.18333D0 + (5143.692D0 / TempK) + (14.613358D0 * log(TempK))
	        A_2   = (21.0894D0 * (Sal ** 0.5D0)) + (0.1248D0 * Sal) - (3.687D-4 * Sal * Sal)
            B_2   = (-772.483D0 * (Sal ** 0.5D0)) - (20.051D0 *Sal)
            C_2   = (-3.3336D0 * (Sal ** 0.5D0))

            !pK2 sigma = 0.011
            pK2   = A_2 + (B_2 / TempK) + (C_2 * log(TempK)) + pK2_0

            K2 = 10.0D0**(-pK2)
        end where
        
   



        where(WhichKs==14);
            ! From Millero, 2010, also for estuarine use.
            ! Marine and Freshwater Research, v. 61, p. 139â€“142.
            ! Fits through compilation of real seawater titration results:
            ! Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
            ! Constants for K's on the SWS;
            ! This is from page 141
            pK10 = -126.34048D0 + (6320.813D0 / TempK) + (19.568224D0 * log(TempK))
            ! This is from their table 2, page 140.
            A1 = (13.4038 * (Sal ** 0.5D0)) + (0.03206D0 * Sal) - (5.242D-5 * Sal *Sal)
            B1 = (-530.659D0 * (Sal ** 0.5D0)) - (5.8210D0 * Sal)
            C1 = (-2.0664D0 * (Sal ** 0.5D0))
            pK1 = pK10 + A1 + (B1 / TempK) + (C1 * log(TempK))  !bug
            K1 = 10.0D0 ** (-pK1)
                 
            ! This is from page 141
            pK20 =  -90.18333D0 + (5143.692D0/TempK) + (14.613358D0 * log(TempK))
            ! This is from their table 3, page 140.
            A2 = (21.3728D0 * (Sal ** 0.5D0)) + (0.1218D0 * Sal) - (3.688D-4 * Sal * Sal)
            B2 = (-788.289D0 * (Sal ** 0.5D0)) - (19.189 * Sal)
            C2 = (-3.374D0 * (Sal ** 0.5D0))
            pK2 = pK20 + A2 + (B2 / TempK) + (C2 * log(TempK))
            K2 = 10.0D0 ** (-pK2)
        end where



        
        !***************************************************************************
        !CorrectKsForPressureNow:
        ! Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
        !       the free scale) are on the SWS scale.
        !       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
        !       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
        !       (the pH scales are the same in this case); the other Ks don't matter.
        !
        !
        ! No salinity dependence is given for the pressure coefficients here.
        ! It is assumed that the salinity is at or very near Sali = 35.
        ! These are valid for the SWS pH scale, but the difference between this and
        ! the total only yields a difference of .004 pH units at 1000 bars, much
        ! less than the uncertainties in the values.
        !****************************************************************************
        ! The sources used are:
        ! Millero, 1995:
        !       Millero, F. J., Thermodynamics of the carbon dioxide system in the
        !       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
        !       See table 9 and eqs. 90-92, p. 675.
        !       TYPO: a factor of 10^3 was left out of the definition of Kappa
        !       TYPO: the value of R given is incorrect with the wrong units
        !       TYPO: the values of the a's for H2S and H2O are from the 1983
        !                values for fresh water
        !       TYPO: the value of a1 for B(OH)3 should be +.1622
        !        Table 9 on p. 675 has no values for Si.
        !       There are a variety of other typos in Table 9 on p. 675.
        !       There are other typos in the paper, and most of the check values
        !       given don't check.
        ! Millero, 1992:
        !       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
        !       CRC Press, 1992. See chapter 6.
        !       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
        !               79, and 96 have typos).
        ! Millero, 1983:
        !       Millero, Frank J., Influence of pressure on chemical processes in
        !       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
        !       Chester, R., Academic Press, 1983.
        !       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
        !       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
        !       these two are necessary to match the values given in Table 43.24
        ! Millero, 1979:
        !       Millero, F. J., The thermodynamics of the carbon dioxide system
        !       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
        !       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
        ! Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
        !       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
        !       This matches the GEOSECS results and is in Edmond and Gieskes.
        ! Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
        !       boric acid, and the pH of seawater, Limnology and Oceanography
        !       13:403-417, 1968.
        ! Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
        !       seawater with respect to calcium carbonate under in situ conditions,
        !       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
        !****************************************************************************
        ! These references often disagree and give different fits for the same thing.
        ! They are not always just an update either; that is, Millero, 1995 may agree
        !       with Millero, 1979, but differ from Millero, 1983.
        ! For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
        !       KP3, and KSi as for the other cases. Peng et al didn't consider the
        !       case of P different from 0. GEOSECS did consider pressure, but didn't
        !       include Phos, Si, or OH, so including the factors here won't matter.
        ! For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
        !       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
        !       including the factors won't matter.
        !****************************************************************************
        !       deltaVs are in cm3/mole
        !       Kappas are in cm3/mole/bar
        !****************************************************************************

        !CorrectK1K2KBForPressure:
        ! MATLAB CODE:
        ! deltaV    = nan(ntps,1); Kappa     = nan(ntps,1);
        ! lnK1fac   = nan(ntps,1); lnK2fac   = nan(ntps,1);
        ! lnKBfac   = nan(ntps,1);
        allocate(deltaV(ntps), Kappa(ntps), lnK1fac(ntps), lnK2fac(ntps), lnKBfac(ntps))

        where (WhichKs==8)
            !***PressureEffectsOnK1inFreshWater:
            !               This is from Millero, 1983.
            deltaV  = -30.54D0 + 0.1849D0 * TempC - 0.0023366D0 * TempC * TempC
            Kappa   = (-6.22D0 + 0.1368D0 * TempC - 0.001233D0  * TempC * TempC) / 1000.0D0;
            lnK1fac = (-deltaV + 0.5D0  * Kappa * Pbar) * Pbar / RT;
            !***PressureEffectsOnK2inFreshWater:
            !               This is from Millero, 1983.
            deltaV  = -29.81D0 + 0.115D0 * TempC - 0.001816D0 * TempC * TempC
            Kappa   = (-5.74D0 + 0.093D0 * TempC - 0.001896D0 * TempC * TempC) / 1000.0D0;
            lnK2fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT
            lnKBfac = 0 ! this doesn't matter since TB = 0 for this case
        end where


        where((WhichKs==6) .or. (WhichKs==7))
            !               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
            !               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
            !               Culberson and Pytkowicz, L and O 13:403-417, 1968:
            !               but the fits are the same as those in
            !               Edmond and Gieskes, GCA, 34:1261-1291, 1970
            !               who in turn quote Li, personal communication
            lnK1fac = (24.2D0 - 0.085D0 * TempC) * Pbar / RT
            lnK2fac = (16.4D0 - 0.04D0  * TempC) * Pbar / RT
            !               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
            !               and matches the GEOSECS results
            lnKBfac = (27.5D0 - 0.095D0 * TempC) * Pbar / RT
        end where


        where ((WhichKs.ne.6) .and. (WhichKs.ne.7) .and. (WhichKs.ne.8))

            !***PressureEffectsOnK1:
            !               These are from Millero, 1995.
            !               They are the same as Millero, 1979 and Millero, 1992.
            !               They are from data of Culberson and Pytkowicz, 1968.

            deltaV  = -25.5D0 + 0.1271D0 * TempC
            !                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979

            Kappa   = (-3.08D0 + 0.0877D0 * TempC) / 1000.0D0
            !                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979

 	        lnK1fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT
            !               The fits given in Millero, 1983 are somewhat different.

            !***PressureEffectsOnK2:
            !               These are from Millero, 1995.
            !               They are the same as Millero, 1979 and Millero, 1992.
            !               They are from data of Culberson and Pytkowicz, 1968.

            deltaV  = -15.82D0 - 0.0219D0 * TempC
            !                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979

            Kappa   = (1.13D0 - 0.1475D0 * TempC) / 1000.0D0
            !                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979

            lnK2fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT
            !               The fit given in Millero, 1983 is different.
            !               Not by a lot for deltaV, but by much for Kappa. %

            !***PressureEffectsOnKB:
            !               This is from Millero, 1979.
            !               It is from data of Culberson and Pytkowicz, 1968.

            deltaV  = -29.48D0 + 0.1622D0 * TempC - 0.002608D0 * TempC * TempC
            !               Millero, 1983 has:
            !                 'deltaV = -28.56 + .1211.*TempCi - .000321.*TempCi.*TempCi
            !              Millero, 1992 has:
            !                 'deltaV = -29.48 + .1622.*TempCi + .295.*(Sali - 34.8)
            !               Millero, 1995 has:
            !                 'deltaV = -29.48 - .1622.*TempCi - .002608.*TempCi.*TempCi
            !                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979

            Kappa = -2.84D0/1000.0D9     ! Millero, 1979
            !               Millero, 1992 and Millero, 1995 also have this.
            !                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
            !               Millero, 1983 has:
            !                 'Kappa = (-3 + .0427.*TempCi)./1000

            lnKBfac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT
        end where

        ! CorrectKWForPressure:
        if(.not.allocated(lnKWfac)) then
            allocate(lnKWfac(ntps))
        endif

        where(WhichKs==8)
            ! PressureEffectsOnKWinFreshWater:
            !               This is from Millero, 1983.
            deltaV  =  -25.6D0 + 0.2324D0 * TempC - 0.0036246D0 * TempC * TempC
            Kappa   = (-7.33D0 + 0.1368D0 * TempC - 0.001233D0  * TempC * TempC) / 1000.0D0
 	        lnKWfac = (-deltaV + 0.5D0 *Kappa * Pbar) * Pbar / RT

            !               NOTE the temperature dependence of KappaK1 and KappaKW
            !               for fresh water in Millero, 1983 are the same.
        end where


        where(WhichKs.ne.8)

            ! GEOSECS doesn't include OH term, so this won't matter.
            ! Peng et al didn't include pressure, but here I assume that the KW correction
            !       is the same as for the other seawater cases.

            ! PressureEffectsOnKW:
            !               This is from Millero, 1983 and his programs CO2ROY(T).BAS.

            deltaV  = -20.02D0 + 0.1119D0 * TempC - 0.001409D0 * TempC * TempC
            !               Millero, 1992 and Millero, 1995 have:

            Kappa   = (-5.13D0 + 0.0794D0 * TempC) / 1000.0D0  !Millero, 1983
            !               Millero, 1995 has this too, but Millero, 1992 is different.

            lnKWfac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT
            !               Millero, 1979 does not list values for these.
        end where

        ! PressureEffectsOnKF:
        !       This is from Millero, 1995, which is the same as Millero, 1983.
        !       It is assumed that KF is on the free pH scale.
        deltaV = -9.78D0 - 0.009D0 * TempC - 0.000942D0 * TempC * TempC
        Kappa = (-3.91D0 + 0.054D0 * TempC) / 1000.0D0
        
        if(.not.(allocated(lnKFfac))) allocate(lnKFfac(ntps))
        lnKFfac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! PressureEffectsOnKS:
        !       This is from Millero, 1995, which is the same as Millero, 1983.
        !       It is assumed that KS is on the free pH scale.
        deltaV = -18.03D0 + 0.0466D0 * TempC + 0.000316D0 * TempC * TempC
        Kappa = (-4.53D0 + 0.09D0 * TempC) / 1000.0D0
        
        if(.not.(allocated(lnKSfac))) allocate(lnKSfac(ntps))
        lnKSfac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! CorrectKP1KP2KP3KSiForPressure:
        ! These corrections don't matter for the GEOSECS choice (WhichKs! = 6) and
        !       the freshwater choice (WhichKs! = 8). For the Peng choice I assume
        !       that they are the same as for the other choices (WhichKs! = 1 to 5).
        ! The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
        !       same as Millero, 1983.

        ! PressureEffectsOnKP1:
        deltaV = -14.51D0 + 0.1211D0 * TempC - 0.000321D0 * TempC * TempC
        Kappa  = (-2.67D0 + 0.0427D0 * TempC) / 1000.0D0

        if(.not.(allocated(lnKP1fac))) allocate(lnKP1fac(ntps))
        lnKP1fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! PressureEffectsOnKP2:
        deltaV = -23.12D0 + 0.1758D0 * TempC - 0.002647D0 * TempC * TempC
        Kappa  = (-5.15D0 + 0.09D0  * TempC) / 1000.0D0

        if(.not.(allocated(lnKP2fac))) allocate(lnKP2fac(ntps))
        lnKP2fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! PressureEffectsOnKP3:
        deltaV = -26.57D0 + 0.202D0 * TempC - 0.003042D0 * TempC * TempC
        Kappa  = (-4.08D0 + 0.0714D0 * TempC) / 1000
        
        if(.not.(allocated(lnKP3fac))) allocate(lnKP3fac(ntps))
        lnKP3fac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! PressureEffectsOnKSi:
        !  The only mention of this is Millero, 1995 where it is stated that the
        !    values have been estimated from the values of boric acid. HOWEVER,
        !    there is no listing of the values in the table.
        !    I used the values for boric acid from above.
        deltaV = -29.48D0 + 0.1622D0 * TempC - 0.002608D0 * TempC * TempC
        Kappa  = -2.84D0 / 1000

        if(.not.(allocated(lnKSifac))) allocate(lnKSifac(ntps))
        lnKSifac = (-deltaV + 0.5D0 * Kappa * Pbar) * Pbar / RT

        ! CorrectKsForPressureHere:
        if(.not.(allocated(K1fac    ))) allocate(K1fac    (ntps))
        if(.not.(allocated(K2fac    ))) allocate(K2fac    (ntps))
        if(.not.(allocated(KWfac    ))) allocate(KWfac    (ntps))
        if(.not.(allocated(KBfac    ))) allocate(KBfac    (ntps))
        if(.not.(allocated(KFfac    ))) allocate(KFfac    (ntps))
        if(.not.(allocated(KSfac    ))) allocate(KSfac    (ntps))
        if(.not.(allocated(KP1fac   ))) allocate(KP1fac   (ntps))
        if(.not.(allocated(KP2fac   ))) allocate(KP2fac   (ntps))
        if(.not.(allocated(KP3fac   ))) allocate(KP3fac   (ntps))
        if(.not.(allocated(KSifac   ))) allocate(KSifac   (ntps))
        if(.not.(allocated(SWStoTOT ))) allocate(SWStoTOT (ntps))
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
        
        K1fac = exp(lnK1fac)
        K1    = K1 * K1fac

        K2fac = exp(lnK2fac)
        K2    = K2 * K2fac
        
        KWfac  = exp(lnKWfac)
        KW     = KW * KWfac
        
        KBfac  = exp(lnKBfac)
        KB     = KB * KBfac
        
        KFfac  = exp(lnKFfac)
        KF     = KF * KFfac
        
        KSfac  = exp(lnKSfac)
        KS     = KS * KSfac
        
        KP1fac = exp(lnKP1fac)
        KP1    = KP1 * KP1fac

        KP2fac = exp(lnKP2fac)
        KP2    = KP2 * KP2fac
        
        KP3fac = exp(lnKP3fac)
        KP3    = KP3 * KP3fac
        
        KSifac = exp(lnKSifac)
        KSi    = KSi * KSifac
        

        ! CorrectpHScaleConversionsForPressure:
        ! fH has been assumed to be independent of pressure.
        SWStoTOT  = (1.0D0 + TS/KS)/(1.0D0 + TS/KS + TF/KF)
        FREEtoTOT =  1.0D0 + TS/KS

        !  The values KS and KF are already pressure-corrected, so the pH scale
        !  conversions are now valid at pressure.

        ! FindpHScaleConversionFactor:
        ! this is the scale they will be put on
        allocate(pHfactor(ntps))

        !Total
        where(pHScale==1)
            pHfactor = SWStoTOT
        end where

        !SWS, they are all on this now
        where(pHScale==2)
            pHfactor = 1.0D0
        end where

        !pHfree
        where(pHScale==3)
            pHfactor = SWStoTOT / FREEtoTOT
        end where

        !pHNBS
        where(pHScale==4)
            pHfactor = fH
        end where

        ! ConvertFromSWSpHScaleToChosenScale:
        K1  = K1  * pHfactor
        K2  = K2  * pHfactor
        KW  = KW  * pHfactor
        KB  = KB  * pHfactor
        KP1 = KP1 * pHfactor
        KP2 = KP2 * pHfactor
        KP3 = KP3 * pHfactor
        KSi = KSi * pHfactor

        ! CalculateFugacityConstants:
        ! This assumes that the pressure is at one atmosphere, or close to it.
        ! Otherwise, the Pres term in the exponent affects the results.
        !       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
        !       Delta and B in cm3/mol
        FugFac = 1.0D0

        if(.not.(allocated(Delta))) allocate(Delta(ntps))
        Delta = 57.7D0 - 0.118D0 * TempK

        if(.not.(allocated(b))) allocate(b(ntps))

        b = -1636.75D0 + (12.0408D0 * TempK) - (0.0327957D0 * TempK * TempK) + &
                (3.16528D0 * 0.00001D0 * TempK * TempK * TempK)

        ! For a mixture of CO2 and air at 1 atm (at low CO2 concentrations)
        P1atm = 1.01325D0 ! in bar
        FugFac = exp((b + 2 * Delta) * P1atm / RT)

        where((WhichKs==6).or.(WhichKs==7))
            FugFac = 1.0D0
        end where

        ! CalculateVPFac:
        ! Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
        !       seawater, Marine Chemistry 8:347-359, 1980.
        ! They fit the data of Goff and Gratch (1946) with the vapor pressure
        !       lowering by sea salt as given by Robinson (1954).
        ! This fits the more complicated Goff and Gratch, and Robinson equations
        !       from 273 to 313 deg K and 0 to 40 Sali with a standard error
        !       of .015!, about 5 uatm over this range.
        ! This may be on IPTS-29 since they didn't mention the temperature scale,
        !       and the data of Goff and Gratch came before IPTS-48.
        ! The references are:
        ! Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
        !       to 212 deg F, Transactions of the American Society of Heating and
        !       Ventilating Engineers 52:95-122, 1946.
        ! Robinson, Journal of the Marine Biological Association of the U. K.
        !       33:449-455, 1954.
        !       This is eq. 10 on p. 350.
        !       This is in atmospheres.
        
        if(.not.(allocated(VPWP))) allocate(VPWP(ntps))
        VPWP = exp(24.4543D0 - 67.4509D0 * (100.0D0 / TempK) - 4.8489D0 * log(TempK / 100.0D0))
        
        if(.not.(allocated(VPCorrWP))) allocate(VPCorrWP(ntps))
        VPCorrWP = exp(-0.000544D0 * Sal)
        
        if(.not.(allocated(VPSWWP))) allocate(VPSWWP(ntps))
        VPSWWP = VPWP * VPCorrWP

        VPFac = 1.0D0 - VPSWWP ! this assumes 1 atmosphere

        !Deallocate  internal arrays. Petras       
                 deallocate(IonS        )     !if(allocated(IonS        ))
                 deallocate(pKS         )     !if(allocated(pKS         ))
                 deallocate(lnKS        )     !if(allocated(lnKS        ))
                 deallocate(lnKF        )     !if(allocated(lnKF        ))
                 deallocate(SWStoTOT    )     !if(allocated(SWStoTOT    ))
                 deallocate(FREEtoTOT   )     !if(allocated(FREEtoTOT   ))
                 deallocate(logKB       )     !if(allocated(logKB       ))
                 deallocate(lnKBtop     )     !if(allocated(lnKBtop     ))
                 deallocate(lnKB        )     !if(allocated(lnKB        ))
                 deallocate(lnKW        )     !if(allocated(lnKW        ))
                 deallocate(lnKP1       )     !if(allocated(lnKP1       ))
                 deallocate(lnKP2       )     !if(allocated(lnKP2       ))
                 deallocate(lnKP3       )     !if(allocated(lnKP3       ))
                 deallocate(lnKSi       )     !if(allocated(lnKSi       ))
                 deallocate(logK1       )     !if(allocated(logK1       ))
                 deallocate(lnK1        )     !if(allocated(lnK1        ))
                 deallocate(pK1         )     !if(allocated(pK1         ))
                 deallocate(logK2       )     !if(allocated(logK2       ))
                 deallocate(lnK2        )     !if(allocated(lnK2        ))
                 deallocate(pK2         )     !if(allocated(pK2         ))
                 deallocate(F1          )     !if(allocated(F1          ))
                 deallocate(F2          )     !if(allocated(F2          ))
                 deallocate(PK1_0       )     !if(allocated(PK1_0       ))
                 deallocate(A_1         )     !if(allocated(A_1         ))
                 deallocate(B_1         )     !if(allocated(B_1         ))
                 deallocate(C_1         )     !if(allocated(C_1         ))
                 deallocate(PK2_0       )     !if(allocated(PK2_0       ))
                 deallocate(A_2         )     !if(allocated(A_2         ))
                 deallocate(B_2         )     !if(allocated(B_2         ))
                 deallocate(C_2         )     !if(allocated(C_2         ))
                 deallocate(PK10        )     !if(allocated(PK10        ))
                 deallocate(A1          )     !if(allocated(A1          ))
                 deallocate(B1          )     !if(allocated(B1          ))
                 deallocate(C1          )     !if(allocated(C1          ))
                 deallocate(PK20        )     !if(allocated(PK20        ))
                 deallocate(A2          )     !if(allocated(A2          ))
                 deallocate(B2          )     !if(allocated(B2          ))
                 deallocate(C2          )     !if(allocated(C2          ))
                 deallocate(deltaV      )     !if(allocated(deltaV      ))
                 deallocate(Kappa       )     !if(allocated(Kappa       ))
                 deallocate(lnK1fac     )     !if(allocated(lnK1fac     ))
                 deallocate(lnK2fac     )     !if(allocated(lnK2fac     ))
                 deallocate(lnKBfac     )     !if(allocated(lnKBfac     ))
                 deallocate(lnKWfac     )     !if(allocated(lnKWfac     ))
                 deallocate(lnKFfac     )     !if(allocated(lnKFfac     ))
                 deallocate(lnKSfac     )     !if(allocated(lnKSfac     ))
                 deallocate(lnKP1fac    )     !if(allocated(lnKP1fac    ))
                 deallocate(lnKP2fac    )     !if(allocated(lnKP2fac    ))
                 deallocate(lnKP3fac    )     !if(allocated(lnKP3fac    ))
                 deallocate(lnKSifac    )     !if(allocated(lnKSifac    ))
                 deallocate(K1fac       )     !if(allocated(K1fac       ))
                 deallocate(K2fac       )     !if(allocated(K2fac       ))
                 deallocate(KWfac       )     !if(allocated(KWfac       ))
                 deallocate(KBfac       )     !if(allocated(KBfac       ))
                 deallocate(KFfac       )     !if(allocated(KFfac       ))
                 deallocate(KSfac       )     !if(allocated(KSfac       ))
                 deallocate(KP1fac      )     !if(allocated(KP1fac      ))
                 deallocate(KP2fac      )     !if(allocated(KP2fac      ))
                 deallocate(KP3fac      )     !if(allocated(KP3fac      ))
                 deallocate(KSifac      )     !if(allocated(KSifac      ))
                 deallocate(pHfactor    )     !if(allocated(pHfactor    ))
                 deallocate(Delta       )     !if(allocated(Delta       ))
                 deallocate(b           )     !if(allocated(b           ))
                 deallocate(VPWP        )     !if(allocated(VPWP        ))
                 deallocate(VPCorrWP    )     !if(allocated(VPCorrWP    ))
                 deallocate(VPSWWP      )     !if(allocated(VPSWWP      ))  
        
        
    end subroutine Constants
    ! **************************************************************************************
    ! End of subroutine that calculates the constants
    ! **************************************************************************************


    subroutine CalculatepHfCO2fromTATC &
               (DEBUG, TAx, TCx, pHx, fCO2x, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)
        ! Outputs pH fCO2, in that order
        ! SUB FindpHfCO2fromTATC, version 01.02, 10-10-97, written by Ernie Lewis.
        ! Inputs: pHScale%, WhichKs%, WhoseKSO4%, TA, TC, Sal, K(), T(), TempC, Pdbar
        ! Outputs: pH, fCO2
        ! This calculates pH and fCO2 from TA and TC at output conditions.

        integer, intent(in) :: DEBUG
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TAx
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: fCO2x

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfCO2fromTATC : Checkpoint #1'
        endif

        ! pH is returned on the scale requested in "pHscale" (see 'constants'...)
        call CalculatepHfromTATC &
             (DEBUG, TAx, TCx, pHx, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfCO2fromTATC : Checkpoint #2'
        endif

        call CalculatefCO2fromTCpH &
             (TCx, pHx, fCO2x, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        !varargout{1} = pHx;
        !varargout{2} = fCO2
        
    end subroutine  CalculatepHfCO2fromTATC


    subroutine CalculatepHfromTATC &
               (DEBUG, TAx, TCx, pHx, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)
        !Outputs pH
        ! SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis.
        ! Inputs: TA, TC, K(), T()
        ! Output: pH
        ! This calculates pH from TA and TC using K1 and K2 by Newton's method.
        ! It tries to solve for the pH at which Residual = 0.
        ! The starting guess is pH = 8.
        ! Though it is coded for H on the total pH scale, for the pH values occuring
        ! in seawater (pH > 6) it will be equally valid on any pH scale (H terms
        ! negligible) as long as the K Constants are on that scale.

        ! Made this to accept vectors. It will continue iterating until all
        ! values in the vector are "abs(deltapH) < pHTol". SVH2007

        integer, intent(in) :: DEBUG
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TAx
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHx

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        ! Auxilary variables?
        real(kind = DBL_PREC), allocatable, dimension(:) :: K1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KWF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3F
        real(kind = DBL_PREC), allocatable, dimension(:) :: TPF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TFF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KFF

        real(kind = DBL_PREC), allocatable, dimension(:) :: deltapH
        integer :: vl
        real(kind = DBL_PREC) :: pHGuess
        real(kind = DBL_PREC) :: pHTol
        real(kind = DBL_PREC) :: ln10
        integer, allocatable, dimension(:) :: INT_VECTOR_1
        integer, allocatable, dimension(:) :: INT_VECTOR_2

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: Denom
        real(kind = DBL_PREC), allocatable, dimension(:) :: CAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: OH
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosTop
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosBot
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfree
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4
        real(kind = DBL_PREC), allocatable, dimension(:) :: HF
        real(kind = DBL_PREC), allocatable, dimension(:) :: Residual
        real(kind = DBL_PREC), allocatable, dimension(:) :: Slope

        integer :: ITERATION_NO
        
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #1'
        endif

        if(.not.(allocated(K1F ))) allocate(K1F (ntps))
        if(.not.(allocated(K2F ))) allocate(K2F (ntps))
        if(.not.(allocated(KWF ))) allocate(KWF (ntps))
        if(.not.(allocated(KP1F))) allocate(KP1F(ntps))
        if(.not.(allocated(KP2F))) allocate(KP2F(ntps))
        if(.not.(allocated(KP3F))) allocate(KP3F(ntps))
        if(.not.(allocated(TPF ))) allocate(TPF (ntps))
        if(.not.(allocated(TSiF))) allocate(TSiF(ntps))
        if(.not.(allocated(KSiF))) allocate(KSiF(ntps))
        if(.not.(allocated(TBF ))) allocate(TBF (ntps))
        if(.not.(allocated(KBF ))) allocate(KBF (ntps))
        if(.not.(allocated(TSF ))) allocate(TSF (ntps))
        if(.not.(allocated(KSF ))) allocate(KSF (ntps))
        if(.not.(allocated(TFF ))) allocate(TFF (ntps))
        if(.not.(allocated(KFF ))) allocate(KFF (ntps))

        K1F  = K1 
        K2F  = K2 
        KWF  = KW 
        KP1F = KP1
        KP2F = KP2
        KP3F = KP3
        TPF  = TP 
        TSiF = TSi
        KSiF = KSi
        TBF  = TB 
        KBF  = KB 
        TSF  = TS 
        KSF  = KS 
        TFF  = TF 
        KFF  = KF 

        if(.not.(allocated(pHx)))     allocate(pHx(ntps))
        if(.not.(allocated(deltapH))) allocate(deltapH(ntps))

        vl      = ntps     ! VectorLength
        pHGuess = 8.0D0    ! this is the first guess
        pHTol   = 0.0001D0 ! tolerance for iterations end
        ln10    = log(10.0D0)

        ! creates a vector holding the first guess for all samples
        ! MATLAB CODE: pHx(1:vl,1) = pHGuess
        pHx     = pHGuess !Program can be imporved by getting pH from last time step
        deltapH = pHTol+1

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #2'
        endif

        !while any(abs(deltapH) > pHTol)
        call LOGICAL_VECTOR_TO_INT_VECTOR(dabs(deltapH) > pHTol, INT_VECTOR_1)

        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3'
        endif

        ITERATION_NO = 0
        
        do while (sum(INT_VECTOR_1) > 0)
            ITERATION_NO = ITERATION_NO + 1
            
            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-1'
            endif

            !MATLAB CODE: H         = 10.^(-pHx);
            if(.not.(allocated(H))) allocate(H(ntps))
            H = 10.0D0 ** (-pHx)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) 'pHx      : ', pHx
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: Denom     = (H.*H + K1F.*H + K1F.*K2F);
            if(.not.(allocated(Denom))) allocate(Denom(ntps))
            Denom  = H * H + K1F * H + K1F * K2F

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Denom    : ', Denom
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) 'K1F      : ', K1F
                write(unit = *, fmt = *) 'K2F      : ', K2F
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: CAlk      = TCx.*K1F.*(H + 2.*K2F)./Denom;
            if(.not.(allocated(CAlk))) allocate(CAlk(ntps))
            CAlk   = TCx * K1F * (H + 2.0D0 * K2F) / Denom

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Calk     : ', Calk
                write(unit = *, fmt = *) 'Denom    : ', Denom
                write(unit = *, fmt = *) 'TCx      : ', TCx
                write(unit = *, fmt = *) 'K1F      : ', K1F
                write(unit = *, fmt = *) 'K2F      : ', K2F
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: BAlk      = TBF.*KBF./(KBF + H);
            if(.not.(allocated(BAlk))) allocate(BAlk(ntps))
            BAlk  = TBF * KBF / (KBF + H)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Balk     : ', Balk
                write(unit = *, fmt = *) 'TBF      : ', TBF
                write(unit = *, fmt = *) 'KBF      : ', KBF
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: OH        = KWF./H;
            if(.not.(allocated(OH))) allocate(OH(ntps))
            OH = KWF / H

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'OH       : ', OH
                write(unit = *, fmt = *) 'KWF      : ', KWF
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
            if(.not.(allocated(PhosTop))) allocate(PhosTop(ntps))
            PhosTop = KP1F * KP2F * H + 2.0D0 * KP1F * KP2F * KP3F - H * H * H

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'PhosTop  : ', PhosTop
                write(unit = *, fmt = *) 'KP1F     : ', KP1F
                write(unit = *, fmt = *) 'KP2F     : ', KP2F
                write(unit = *, fmt = *) 'KP3F     : ', KP3F
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
            if(.not.(allocated(PhosBot))) allocate(PhosBot(ntps))
            PhosBot = H * H * H + KP1F * H * H + KP1F * KP2F * H + KP1F * KP2F * KP3F

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'PhosTop  : ', PhosBot
                write(unit = *, fmt = *) 'KP1F     : ', KP1F
                write(unit = *, fmt = *) 'KP2F     : ', KP2F
                write(unit = *, fmt = *) 'KP3F     : ', KP3F
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: PAlk      = TPF.*PhosTop./PhosBot;
            if(.not.(allocated(PAlk))) allocate(PAlk(ntps))
            PAlk   = TPF * PhosTop / PhosBot

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'PAlk     : ', PAlk
                write(unit = *, fmt = *) 'TPF      : ', TPF
                write(unit = *, fmt = *) 'PhosTop  : ', PhosTop
                write(unit = *, fmt = *) 'PhosTop  : ', PhosBot
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: SiAlk     = TSiF.*KSiF./(KSiF + H)
            if(.not.(allocated(SiAlk))) allocate(SiAlk(ntps))
            SiAlk  = TSiF * KSiF / (KSiF + H)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'SiAlk    : ', SiAlk
                write(unit = *, fmt = *) 'TSiF     : ', TSiF
                write(unit = *, fmt = *) 'KSiF     : ', TSiF
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            ! pH scale conversion factor
            !MATLAB CODE: FREEtoTOT = (1 + TSF./KSF); % pH scale conversion factor
            if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
            FREEtoTOT = 1.0D0 + TSF / KSF

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'FREEtoTOT: ', FREEtoTOT
                write(unit = *, fmt = *) 'TSF      : ', TSF
                write(unit = *, fmt = *) 'KSF      : ', KSF
                write(unit = *, fmt = *) '--------------------'
            endif

            ! for H on the total scale
            !MATLAB CODE: Hfree     = H./FREEtoTOT; % for H on the total scale
            if(.not.(allocated(Hfree))) allocate(Hfree(ntps))
            call ASSIGN_DBL_VECTOR_CONTENT(Hfree    , H / FREEtoTOT)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Hfree    : ', Hfree
                write(unit = *, fmt = *) 'FREEtoTOT: ', FREEtoTOT
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) '--------------------'
            endif

            ! since KS is on the free scale
            !MATLAB CODE: HSO4      = TSF./(1 + KSF./Hfree); % since KS is on the free scale
            if(.not.(allocated(HSO4))) allocate(HSO4(ntps))
            HSO4     = TSF / (1.0D0 + KSF / Hfree)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'HSO4     : ', HSO4
                write(unit = *, fmt = *) 'TSF      : ', TSF
                write(unit = *, fmt = *) 'KSF      : ', KSF
                write(unit = *, fmt = *) 'Hfree    : ', Hfree
                write(unit = *, fmt = *) '--------------------'
            endif

            ! since KF is on the free scale
            !MATLAB CODE: HF        = TFF./(1 + KFF./Hfree); % since KF is on the free scale
            if(.not.(allocated(HF))) allocate(HF(ntps))
            HF       = TFF / (1.0D0 + KFF / Hfree)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'HF       : ', HF
                write(unit = *, fmt = *) 'KFF      : ', KFF
                write(unit = *, fmt = *) 'Hfree    : ', Hfree
                write(unit = *, fmt = *) '--------------------'
            endif

            !MATLAB CODE: Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF;
            if(.not.(allocated(Residual))) allocate(Residual(ntps))
            Residual = TAx - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Residual : ', Residual
                write(unit = *, fmt = *) 'TAx      : ', TAx
                write(unit = *, fmt = *) 'CAlk     : ', CAlk
                write(unit = *, fmt = *) 'BAlk     : ', BAlk
                write(unit = *, fmt = *) 'OH       : ', OH
                write(unit = *, fmt = *) 'PAlk     : ', PAlk
                write(unit = *, fmt = *) 'SiAlk    : ', SiAlk
                write(unit = *, fmt = *) 'Hfree    : ', Hfree
                write(unit = *, fmt = *) 'HSO4     : ', HSO4
                write(unit = *, fmt = *) 'HF       : ', HF
                write(unit = *, fmt = *) '--------------------'
            endif

            ! find Slope dTA/dpH;
            ! (this is not exact, but keeps all important terms);
            !MATLAB CODE: Slope     = ln10.*(TCx.*K1F.*H.*(H.*H + K1F.*K2F + 4.*H.*K2F)./Denom./Denom + BAlk.*H./(KBF + H) + OH + H);
            if(.not.(allocated(Slope))) allocate(Slope(ntps))
            
            Slope    = ln10 * &
                 (TCx * K1F * H * (H * H + K1F * K2F + 4.0D0 * H * K2F) / Denom / Denom + BAlk * H / (KBF + H) + OH + H)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'Slope    : ', Slope
                write(unit = *, fmt = *) 'TCx      : ', TCx
                write(unit = *, fmt = *) 'CAlk     : ', CAlk
                write(unit = *, fmt = *) 'K1F      : ', K1F
                write(unit = *, fmt = *) 'K2F      : ', K2F
                write(unit = *, fmt = *) 'Denom    : ', Denom
                write(unit = *, fmt = *) 'BAlk     : ', BAlk
                write(unit = *, fmt = *) 'KBF      : ', KBF
                write(unit = *, fmt = *) 'H        : ', H
                write(unit = *, fmt = *) 'OH       : ', OH
                write(unit = *, fmt = *) '--------------------'
            endif

            ! this is Newton's method
            !MATLAB CODE: deltapH   = Residual./Slope; % this is Newton's method
            if(.not.(allocated(deltapH))) allocate(deltapH(ntps))
            deltapH  = Residual / Slope

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'deltapH  : ', deltapH
                write(unit = *, fmt = *) 'Residual : ', Residual
                write(unit = *, fmt = *) 'Slope    : ', Slope
                write(unit = *, fmt = *) '********************'
                write(unit = *, fmt = *) '********************'
                write(unit = *, fmt = *) '********************'
            endif

            ! to keep the jump from being too big
            !MATLAB CODE : while any(abs(deltapH) > 1)
            call LOGICAL_VECTOR_TO_INT_VECTOR(dabs(deltapH) > 1.0D0, INT_VECTOR_2)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-2'
            endif

            do while (sum(INT_VECTOR_2) > 0)
                !MATLAB CODE:
                !FF=abs(deltapH)>1
                !deltapH(FF)=deltapH(FF)./2;
                if (DEBUG .eq. 1) then
                    write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-2-1'
                    write(unit = *, fmt = *) 'deltapH : ', deltapH
                endif

                where(dabs(deltapH) > 1)
                    deltapH = deltapH / 2.0D0
                end where

                if (DEBUG .eq. 1) then
                    write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-2-2'
                    write(unit = *, fmt = *) 'deltapH : ', deltapH
                endif

                call LOGICAL_VECTOR_TO_INT_VECTOR(dabs(deltapH) > 1, INT_VECTOR_2)

            end do

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-3'
                write(unit = *, fmt = *) INT_VECTOR_2
                write(*,*) 'deltapH : ', deltapH
            endif
            ! Is on the same scale as K1 and K2 were calculated...
            pHx = pHx + deltapH

            call LOGICAL_VECTOR_TO_INT_VECTOR(dabs(deltapH) > pHTol, INT_VECTOR_1)

            if (DEBUG .eq. 1) then
                write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #3-4'
                write(*,*) INT_VECTOR_1
                write(*,*) 'deltapH : ', deltapH
                write(*,*) 'pHTol   : ', pHTol
            end if

            if (ITERATION_NO > 1000) then
                write(*,*) 'PH DUMP'
                write(*,*) pHx
                write(*,*) 'deltapH DUMP'
                write(*,*) deltapH
                stop 
            end if
        end do

        !write(*,*)  'NUMBER OF ITERATIONS NEEDED = ', ITERATION_NO
        
        if (DEBUG .eq. 1) then
            write(unit = *, fmt = *) 'subroutine CalculatepHfromTATC : Checkpoint #4'
        endif
        !varargout{1}=pHx;
        
        
!         Deallocate auxilary arrays. Petras

        deallocate (K1F            )     !if(allocated(K1F          ))
        deallocate (K2F            )     !if(allocated(K2F          ))
        deallocate (KWF            )     !if(allocated(KWF          ))
        deallocate (KP1F           )     !if(allocated(KP1F         ))
        deallocate (KP2F           )     !if(allocated(KP2F         ))
        deallocate (KP3F           )     !if(allocated(KP3F         ))
        deallocate (TPF            )     !if(allocated(TPF          ))
        deallocate (TSiF           )     !if(allocated(TSiF         ))
        deallocate (KSiF           )     !if(allocated(KSiF         ))
        deallocate (TBF            )     !if(allocated(TBF          ))
        deallocate (KBF            )     !if(allocated(KBF          ))
        deallocate (TSF            )     !if(allocated(TSF          ))
        deallocate (KSF            )     !if(allocated(KSF          ))
        deallocate (TFF            )     !if(allocated(TFF          ))
        deallocate (KFF            )     !if(allocated(KFF          ))
        deallocate (deltapH        )     !if(allocated(deltapH      ))
        deallocate (INT_VECTOR_1   )     !if(allocated(INT_VECTOR_1 ))
        deallocate (INT_VECTOR_2   )     !if(allocated(INT_VECTOR_2 ))
        deallocate (H              )     !if(allocated(H            ))
        deallocate (Denom          )     !if(allocated(Denom        ))
        deallocate (CAlk           )     !if(allocated(CAlk         ))
        deallocate (BAlk           )     !if(allocated(BAlk         ))
        deallocate (OH             )     !if(allocated(OH           ))
        deallocate (PhosTop        )     !if(allocated(PhosTop      ))
        deallocate (PhosBot        )     !if(allocated(PhosBot      ))
        deallocate (PAlk           )     !if(allocated(PAlk         ))
        deallocate (SiAlk          )     !if(allocated(SiAlk        ))
        deallocate (FREEtoTOT      )     !if(allocated(FREEtoTOT    ))
        deallocate (Hfree          )     !if(allocated(Hfree        ))
        deallocate (HSO4           )     !if(allocated(HSO4         ))
        deallocate (HF             )     !if(allocated(HF           ))
        deallocate (Residual       )     !if(allocated(Residual     ))
        deallocate (Slope          )     !if(allocated(Slope        ))
        
        
    end subroutine CalculatepHfromTATC


    subroutine CalculatefCO2fromTCpH &
               (TCx, pHx, fCO2x, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculatefCO2fromTCpH, version 02.02, 12-13-96, written by Ernie Lewis.
        ! ' Inputs: TC, pH, K0, K1, K2
        ! ' Output: fCO2
        ! ' This calculates fCO2 from TC and pH, using K0, K1, and K2.

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCx
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pHx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: fCO2x

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS


        real(kind = DBL_PREC), allocatable, dimension(:) :: H

        if(.not.(allocated(H))) allocate(H(ntps))
        H     = 10.0D0 ** (-pHx)

        if(.not.(allocated(fCO2x))) allocate(fCO2x(ntps))
        fCO2x = TCx * H * H / (H * H + K1 * H + K1 * K2 ) / K0
        !varargout{1} = fCO2x;
        
        deallocate(H)
        
    end subroutine CalculatefCO2fromTCpH



    subroutine CalculateTCfromTApH &
               (TAx, pHx, TCxtemp, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculateTCfromTApH, version 02.03, 10-10-97, written by Ernie Lewis.
        ! ' Inputs: TA, pH, K(), T()
        ! ' Output: TC
        ! ' This calculates TC from TA and pH.
        ! ' Though it is coded for H on the total pH scale, for the pH values occuring
        ! ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
        ! ' negligible) as long as the K Constants are on that scale.

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TAx
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pHx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: TCxtemp

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: K1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KWF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3F
        real(kind = DBL_PREC), allocatable, dimension(:) :: TPF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TFF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KFF

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: OH
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosTop
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosBot
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfree
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4
        real(kind = DBL_PREC), allocatable, dimension(:) :: HF
        real(kind = DBL_PREC), allocatable, dimension(:) :: CAlk

        if(.not.(allocated(K1F ))) allocate(K1F (ntps))
        if(.not.(allocated(K2F ))) allocate(K2F (ntps))
        if(.not.(allocated(KWF ))) allocate(KWF (ntps))
        if(.not.(allocated(KP1F))) allocate(KP1F(ntps))
        if(.not.(allocated(KP2F))) allocate(KP2F(ntps))
        if(.not.(allocated(KP3F))) allocate(KP3F(ntps))
        if(.not.(allocated(TPF ))) allocate(TPF (ntps))
        if(.not.(allocated(TSiF))) allocate(TSiF(ntps))
        if(.not.(allocated(KSiF))) allocate(KSiF(ntps))
        if(.not.(allocated(TBF ))) allocate(TBF (ntps))
        if(.not.(allocated(KBF ))) allocate(KBF (ntps))
        if(.not.(allocated(TSF ))) allocate(TSF (ntps))
        if(.not.(allocated(KSF ))) allocate(KSF (ntps))
        if(.not.(allocated(TFF ))) allocate(TFF (ntps))
        if(.not.(allocated(KFF ))) allocate(KFF (ntps))
        
        K1F  = K1
        K2F  = K2
        KWF  = KW
        KP1F = KP1
        KP2F = KP2
        KP3F = KP3
        TPF  = TP
        TSiF = TSi
        KSiF = KSi
        TBF  = TB
        KBF  = KB
        TSF  = TS
        KSF  = KS
        TFF  = TF
        KFF  = KF

        
        if(.not.(allocated(H      ))) allocate(H      (ntps))
        if(.not.(allocated(BAlk   ))) allocate(BAlk   (ntps))
        if(.not.(allocated(OH     ))) allocate(OH     (ntps))
        if(.not.(allocated(PhosTop))) allocate(PhosTop(ntps))
        if(.not.(allocated(PhosBot))) allocate(PhosBot(ntps))
        if(.not.(allocated(PAlk   ))) allocate(PAlk   (ntps))
        if(.not.(allocated(SiAlk  ))) allocate(SiAlk  (ntps))
        
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
        if(.not.(allocated(Hfree    ))) allocate(Hfree    (ntps))
        if(.not.(allocated(HSO4     ))) allocate(HSO4     (ntps))
        if(.not.(allocated(HF       ))) allocate(HF       (ntps))
        if(.not.(allocated(CAlk     ))) allocate(CAlk     (ntps))
        if(.not.(allocated(TCxtemp  ))) allocate(TCxtemp  (ntps))
                
        H        = 10.0D0 ** (-pHx)
        BAlk     = TBF * KBF / (KBF + H)
        OH       = KWF / H
        PhosTop  = KP1F * KP2F * H + 2.0D0 * KP1F * KP2F * KP3F - H * H * H
        PhosBot  = H * H * H + KP1F * H * H + KP1F * KP2F * H + KP1F * KP2F * KP3F
        PAlk     = TPF * PhosTop / PhosBot
        SiAlk    = TSiF * KSiF / (KSiF + H)

        ! pH scale conversion factor
        FREEtoTOT = (1.0D0 + TSF / KSF)

        !' for H on the total scale
        Hfree    = H / FREEtoTOT

        !' since KS is on the free scale
        HSO4     = TSF / (1.0D0 + KSF / Hfree)

        !' since KF is on the free scale
        HF       = TFF / (1.0D0 + KFF / Hfree)
        CAlk     = TAx - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF
        TCxtemp  = CAlk * (H * H + K1F * H + K1F * K2F) / (K1F * (H + 2.0D0 * K2F))
        !varargout{1} = TCxtemp;
        
!        Deallocate auxilary arrays. Petras        
          deallocate(K1F      )    !if(allocated(K1F      ))
          deallocate(K2F      )    !if(allocated(K2F      ))
          deallocate(KWF      )    !if(allocated(KWF      ))
          deallocate(KP1F     )    !if(allocated(KP1F     ))
          deallocate(KP2F     )    !if(allocated(KP2F     ))
          deallocate(KP3F     )    !if(allocated(KP3F     ))
          deallocate(TPF      )    !if(allocated(TPF      ))
          deallocate(TSiF     )    !if(allocated(TSiF     ))
          deallocate(KSiF     )    !if(allocated(KSiF     ))
          deallocate(TBF      )    !if(allocated(TBF      ))
          deallocate(KBF      )    !if(allocated(KBF      ))
          deallocate(TSF      )    !if(allocated(TSF      ))
          deallocate(KSF      )    !if(allocated(KSF      ))
          deallocate(TFF      )    !if(allocated(TFF      ))
          deallocate(KFF      )    !if(allocated(KFF      ))
          deallocate(H        )    !if(allocated(H        ))
          deallocate(BAlk     )    !if(allocated(BAlk     ))
          deallocate(OH       )    !if(allocated(OH       ))
          deallocate(PhosTop  )    !if(allocated(PhosTop  ))
          deallocate(PhosBot  )    !if(allocated(PhosBot  ))
          deallocate(PAlk     )    !if(allocated(PAlk     ))
          deallocate(SiAlk    )    !if(allocated(SiAlk    ))
          deallocate(FREEtoTOT)    !if(allocated(FREEtoTOT))
          deallocate(Hfree    )    !if(allocated(Hfree    ))
          deallocate(HSO4     )    !if(allocated(HSO4     ))
          deallocate(HF       )    !if(allocated(HF       ))
          deallocate(CAlk     )    !if(allocated(CAlk     ))
          
    end subroutine CalculateTCfromTApH


    subroutine CalculatepHfromTAfCO2 &
               (TAi, fCO2i, pH, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculatepHfromTAfCO2, version 04.01, 10-13-97, written by Ernie Lewis.
        ! ' Inputs: TA, fCO2, K0, K(), T()
        ! ' Output: pH
        ! ' This calculates pH from TA and fCO2 using K1 and K2 by Newton's method.
        ! ' It tries to solve for the pH at which Residual = 0.
        ! ' The starting guess is pH = 8.
        ! ' Though it is coded for H on the total pH scale, for the pH values occuring
        ! ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
        ! ' negligible) as long as the K Constants are on that scale.

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TAi
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: fCO2i
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pH

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: K0F
        real(kind = DBL_PREC), allocatable, dimension(:) :: K1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KWF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3F
        real(kind = DBL_PREC), allocatable, dimension(:) :: TPF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TFF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KFF

        real(kind = DBL_PREC), allocatable, dimension(:) :: deltapH
        integer :: vl
        real(kind = DBL_PREC) :: pHGuess
        real(kind = DBL_PREC) :: pHTol
        real(kind = DBL_PREC) :: ln10
        integer, allocatable, dimension(:) :: INT_VECTOR_1
        integer, allocatable, dimension(:) :: INT_VECTOR_2

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: HCO3
        real(kind = DBL_PREC), allocatable, dimension(:) :: CO3
        real(kind = DBL_PREC), allocatable, dimension(:) :: Denom
        real(kind = DBL_PREC), allocatable, dimension(:) :: CAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: OH
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosTop
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosBot
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfree
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4
        real(kind = DBL_PREC), allocatable, dimension(:) :: HF
        real(kind = DBL_PREC), allocatable, dimension(:) :: Residual
        real(kind = DBL_PREC), allocatable, dimension(:) :: Slope

        
        if(.not.(allocated(K0F )))  allocate(K0F (ntps))
        if(.not.(allocated(K1F )))  allocate(K1F (ntps))
        if(.not.(allocated(K2F )))  allocate(K2F (ntps))
        if(.not.(allocated(KWF )))  allocate(KWF (ntps))
        if(.not.(allocated(KP1F)))  allocate(KP1F(ntps))
        if(.not.(allocated(KP2F)))  allocate(KP2F(ntps))
        if(.not.(allocated(KP3F)))  allocate(KP3F(ntps))
        if(.not.(allocated(TPF )))  allocate(TPF (ntps))
        if(.not.(allocated(TSiF)))  allocate(TSiF(ntps))
        if(.not.(allocated(KSiF)))  allocate(KSiF(ntps))
        if(.not.(allocated(TBF )))  allocate(TBF (ntps))
        if(.not.(allocated(KBF )))  allocate(KBF (ntps))
        if(.not.(allocated(TSF )))  allocate(TSF (ntps))
        if(.not.(allocated(KSF )))  allocate(KSF (ntps))
        if(.not.(allocated(TFF )))  allocate(TFF (ntps))
        if(.not.(allocated(KFF )))  allocate(KFF (ntps))

        K0F  = K0 
        K1F  = K1 
        K2F  = K2 
        KWF  = KW 
        KP1F = KP1
        KP2F = KP2
        KP3F = KP3
        TPF  = TP 
        TSiF = TSi
        KSiF = KSi
        TBF  = TB 
        KBF  = KB 
        TSF  = TS 
        KSF  = KS 
        TFF  = TF 
        KFF  = KF 

        vl          = ntps    ! VectorLength
        pHGuess     = 8.0D0   ! this is the first guess
        pHTol       = 0.0001  ! tolerance for iterations end
        ln10        = log(10.0D0)

        ! creates a vector holding the first guess for all samples
        ! MATLAB CODE: pHx(1:vl,1) = pHGuess
        pH      = pHGuess !Program can be imporved by getting pH from last time step
        deltapH = pHTol+1


        !while any(abs(deltapH) > pHTol)
        call LOGICAL_VECTOR_TO_INT_VECTOR(abs(deltapH) > pHTol, INT_VECTOR_1)

        do while (sum(INT_VECTOR_1) > 0)
            if(.not.(allocated(H        ))) allocate(H        (ntps))
            if(.not.(allocated(HCO3     ))) allocate(HCO3     (ntps))
            if(.not.(allocated(CO3      ))) allocate(CO3      (ntps))
            if(.not.(allocated(CAlk     ))) allocate(CAlk     (ntps))
            if(.not.(allocated(BAlk     ))) allocate(BAlk     (ntps))
            if(.not.(allocated(OH       ))) allocate(OH       (ntps))
            if(.not.(allocated(PhosTop  ))) allocate(PhosTop  (ntps))
            if(.not.(allocated(PhosBot  ))) allocate(PhosBot  (ntps))
            if(.not.(allocated(PAlk     ))) allocate(PAlk     (ntps))
            if(.not.(allocated(SiAlk    ))) allocate(SiAlk    (ntps))
            if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
            if(.not.(allocated(Hfree    ))) allocate(Hfree    (ntps))
            if(.not.(allocated(HSO4     ))) allocate(HSO4     (ntps))
            if(.not.(allocated(HF       ))) allocate(HF       (ntps))
            if(.not.(allocated(Residual ))) allocate(Residual (ntps))
            if(.not.(allocated(Slope    ))) allocate(Slope    (ntps))
            if(.not.(allocated(deltapH  ))) allocate(deltapH  (ntps))

            call ASSIGN_DBL_VECTOR_CONTENT(H        , 10.D0 ** (-pH))
            call ASSIGN_DBL_VECTOR_CONTENT(HCO3     , K0F * K1F * fCO2i / H)
            call ASSIGN_DBL_VECTOR_CONTENT(CO3      , K0F * K1F * K2F * fCO2i / (H * H))
            call ASSIGN_DBL_VECTOR_CONTENT(CAlk     , HCO3 + 2.0D0 * CO3)
            call ASSIGN_DBL_VECTOR_CONTENT(BAlk     , TBF * KBF / (KBF + H))
            call ASSIGN_DBL_VECTOR_CONTENT(OH       , KWF / H)
            call ASSIGN_DBL_VECTOR_CONTENT(PhosTop  , KP1F * KP2F * H + 2.0D0 * KP1F * KP2F * KP3F - H * H *H)
            call ASSIGN_DBL_VECTOR_CONTENT(PhosBot  , H * H * H + KP1F * H * H + KP1F * KP2F * H + KP1F * KP2F * KP3F)
            call ASSIGN_DBL_VECTOR_CONTENT(PAlk     , TPF * PhosTop / PhosBot)
            call ASSIGN_DBL_VECTOR_CONTENT(SiAlk    , TSiF * KSiF / (KSiF + H))

            ! ' pH scale conversion factor
            call ASSIGN_DBL_VECTOR_CONTENT(FREEtoTOT, (1.0D0 + TSF / KSF))

            !' for H on the total scale
            call ASSIGN_DBL_VECTOR_CONTENT(Hfree    , H / FREEtoTOT)

            !' since KS is on the free scale
            call ASSIGN_DBL_VECTOR_CONTENT(HSO4     , TSF / (1.0D0 + KSF/Hfree))

            ! ' since KF is on the free scale
            call ASSIGN_DBL_VECTOR_CONTENT(HF       , TFF / (1.0D0 + KFF/Hfree))
            call ASSIGN_DBL_VECTOR_CONTENT(Residual , TAi - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF)

            ! find Slope dTA/dpH (this is not exact, but keeps all important terms):
            call ASSIGN_DBL_VECTOR_CONTENT(Slope    , ln10 * (HCO3 + 4.0D0 * CO3 + BAlk * H / (KBF + H) + OH + H))
            call ASSIGN_DBL_VECTOR_CONTENT(deltapH  , Residual / Slope)

            ! ' to keep the jump from being too big:
            !MATLAB CODE: while any(abs(deltapH) > 1)
            call LOGICAL_VECTOR_TO_INT_VECTOR(abs(deltapH) > 1, INT_VECTOR_2)

            do while (sum(INT_VECTOR_2) > 0)
                !MATLAB CODE:
                !FF=abs(deltapH)>1
                !deltapH(FF)=deltapH(FF)./2;
                where(abs(deltapH) > 1)
                    deltapH = deltapH / 2
                end where

                call LOGICAL_VECTOR_TO_INT_VECTOR(abs(deltapH) > 1, INT_VECTOR_2)
            end do

            pH = pH + deltapH

            call LOGICAL_VECTOR_TO_INT_VECTOR(abs(deltapH) > pHTol, INT_VECTOR_1)
        end do

    !varargout{1}=pH;
    
!     Deallocate auxilary arrayas. Petras    
        deallocate(K0F         )    !if(allocated(K0F           ))
        deallocate(K1F         )    !if(allocated(K1F           ))
        deallocate(K2F         )    !if(allocated(K2F           ))
        deallocate(KWF         )    !if(allocated(KWF           ))
        deallocate(KP1F        )    !if(allocated(KP1F          ))
        deallocate(KP2F        )    !if(allocated(KP2F          ))
        deallocate(KP3F        )    !if(allocated(KP3F          ))
        deallocate(TPF         )    !if(allocated(TPF           ))
        deallocate(TSiF        )    !if(allocated(TSiF          ))
        deallocate(KSiF        )    !if(allocated(KSiF          ))
        deallocate(TBF         )    !if(allocated(TBF           ))
        deallocate(KBF         )    !if(allocated(KBF           ))
        deallocate(TSF         )    !if(allocated(TSF           ))
        deallocate(KSF         )    !if(allocated(KSF           ))
        deallocate(TFF         )    !if(allocated(TFF           ))
        deallocate(KFF         )    !if(allocated(KFF           ))
        deallocate(deltapH     )    !if(allocated(deltapH       ))
        deallocate(INT_VECTOR_1)    !if(allocated(INT_VECTOR_1  ))
        deallocate(INT_VECTOR_2)    !if(allocated(INT_VECTOR_2  ))
        deallocate(H           )    !if(allocated(H             ))
        deallocate(HCO3        )    !if(allocated(HCO3          ))
        deallocate(CO3         )    !if(allocated(CO3           ))
        deallocate(Denom       )    !if(allocated(Denom         ))
        deallocate(CAlk        )    !if(allocated(CAlk          ))
        deallocate(BAlk        )    !if(allocated(BAlk          ))
        deallocate(OH          )    !if(allocated(OH            ))
        deallocate(PhosTop     )    !if(allocated(PhosTop       ))
        deallocate(PhosBot     )    !if(allocated(PhosBot       ))
        deallocate(PAlk        )    !if(allocated(PAlk          ))
        deallocate(SiAlk       )    !if(allocated(SiAlk         ))
        deallocate(FREEtoTOT   )    !if(allocated(FREEtoTOT     ))
        deallocate(Hfree       )    !if(allocated(Hfree         ))
        deallocate(HSO4        )    !if(allocated(HSO4          ))
        deallocate(HF          )    !if(allocated(HF            ))
        deallocate(Residual    )    !if(allocated(Residual      ))
        deallocate(Slope       )    !if(allocated(Slope         ))
        
    end subroutine CalculatepHfromTAfCO2


    subroutine CalculateTAfromTCpH &
               (TCi, pHi, TActemp, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculateTAfromTCpH, version 02.02, 10-10-97, written by Ernie Lewis.
        ! ' Inputs: TC, pH, K(), T()
        ! ' Output: TA
        ! ' This calculates TA from TC and pH.
        ! ' Though it is coded for H on the total pH scale, for the pH values occuring
        ! ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
        ! ' negligible) as long as the K Constants are on that scale.

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCi
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pHi
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: TActemp

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: K1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: K2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KWF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP1F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP2F
        real(kind = DBL_PREC), allocatable, dimension(:) :: KP3F
        real(kind = DBL_PREC), allocatable, dimension(:) :: TPF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSiF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KBF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KSF
        real(kind = DBL_PREC), allocatable, dimension(:) :: TFF
        real(kind = DBL_PREC), allocatable, dimension(:) :: KFF

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: BAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: OH
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosTop
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosBot
        real(kind = DBL_PREC), allocatable, dimension(:) :: PAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: SiAlk
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: Hfree
        real(kind = DBL_PREC), allocatable, dimension(:) :: HSO4
        real(kind = DBL_PREC), allocatable, dimension(:) :: HF
        real(kind = DBL_PREC), allocatable, dimension(:) :: CAlk

        
        if(.not.(allocated(K1F )))  allocate(K1F (ntps))
        if(.not.(allocated(K2F )))  allocate(K2F (ntps))
        if(.not.(allocated(KWF )))  allocate(KWF (ntps))
        if(.not.(allocated(KP1F)))  allocate(KP1F(ntps))
        if(.not.(allocated(KP2F)))  allocate(KP2F(ntps))
        if(.not.(allocated(KP3F)))  allocate(KP3F(ntps))
        if(.not.(allocated(TPF )))  allocate(TPF (ntps))
        if(.not.(allocated(TSiF)))  allocate(TSiF(ntps))
        if(.not.(allocated(KSiF)))  allocate(KSiF(ntps))
        if(.not.(allocated(TBF )))  allocate(TBF (ntps))
        if(.not.(allocated(KBF )))  allocate(KBF (ntps))
        if(.not.(allocated(TSF )))  allocate(TSF (ntps))
        if(.not.(allocated(KSF )))  allocate(KSF (ntps))
        if(.not.(allocated(TFF )))  allocate(TFF (ntps))
        if(.not.(allocated(KFF )))  allocate(KFF (ntps))
        
        K1F  = K1
        K2F  = K2
        KWF  = KW
        KP1F = KP1
        KP2F = KP2
        KP3F = KP3
        TPF  = TP
        TSiF = TSi
        KSiF = KSi
        TBF  = TB
        KBF  = KB
        TSF  = TS
        KSF  = KS
        TFF  = TF
        KFF  = KF

        if(.not.(allocated(H        ))) allocate(H        (ntps))
        if(.not.(allocated(CAlk     ))) allocate(CAlk     (ntps))
        if(.not.(allocated(BAlk     ))) allocate(BAlk     (ntps))
        if(.not.(allocated(OH       ))) allocate(OH       (ntps))
        if(.not.(allocated(PhosTop  ))) allocate(PhosTop  (ntps))
        if(.not.(allocated(PhosBot  ))) allocate(PhosBot  (ntps))
        if(.not.(allocated(PAlk     ))) allocate(PAlk     (ntps))
        if(.not.(allocated(SiAlk    ))) allocate(SiAlk    (ntps))
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
        if(.not.(allocated(Hfree    ))) allocate(Hfree    (ntps))
        if(.not.(allocated(HSO4     ))) allocate(HSO4     (ntps))
        if(.not.(allocated(HF       ))) allocate(HF       (ntps))
        if(.not.(allocated(TActemp  ))) allocate(TActemp  (ntps))
        
        H        = 10.0D0 ** (-pHi)
        CAlk     = TCi * K1F * (H + 2.0D0 * K2F) / (H * H + K1F * H + K1F * K2F)
        BAlk     = TBF * KBF / (KBF + H)
        OH       = KWF / H
        PhosTop  = KP1F * KP2F * H + 2.0D0 * KP1F * KP2F * KP3F - H * H * H
        PhosBot  = H * H * H + KP1F * H * H + KP1F * KP2F * H + KP1F * KP2F * KP3F
        PAlk     = TPF * PhosTop / PhosBot
        SiAlk    = TSiF * KSiF / (KSiF + H)

        ! ' pH scale conversion factor
        FREEtoTOT = (1.0D0 + TSF / KSF)

        !' for H on the total scale
        Hfree    = H / FREEtoTOT

        ! ' since KS is on the free scale
        HSO4     = TSF /(1.0D0 + KSF / Hfree)

        ! ' since KF is on the free scale
        HF       = TFF / (1.0D0 + KFF / Hfree)
        TActemp  = CAlk + BAlk + OH + PAlk + SiAlk - Hfree - HSO4 - HF
        !varargout{1}=TActemp;

        !Deallocate auxillary arrays. Petras               
          deallocate(K1F        )   !if(allocated(K1F      ))
          deallocate(K2F        )   !if(allocated(K2F      ))
          deallocate(KWF        )   !if(allocated(KWF      ))
          deallocate(KP1F       )   !if(allocated(KP1F     ))
          deallocate(KP2F       )   !if(allocated(KP2F     ))
          deallocate(KP3F       )   !if(allocated(KP3F     ))
          deallocate(TPF        )   !if(allocated(TPF      ))
          deallocate(TSiF       )   !if(allocated(TSiF     ))
          deallocate(KSiF       )   !if(allocated(KSiF     ))
          deallocate(TBF        )   !if(allocated(TBF      ))
          deallocate(KBF        )   !if(allocated(KBF      ))
          deallocate(TSF        )   !if(allocated(TSF      ))
          deallocate(KSF        )   !if(allocated(KSF      ))
          deallocate(TFF        )   !if(allocated(TFF      ))
          deallocate(KFF        )   !if(allocated(KFF      ))
          deallocate(H          )   !if(allocated(H        ))
          deallocate(CAlk       )   !if(allocated(CAlk     ))
          deallocate(BAlk       )   !if(allocated(BAlk     ))
          deallocate(OH         )   !if(allocated(OH       ))
          deallocate(PhosTop    )   !if(allocated(PhosTop  ))
          deallocate(PhosBot    )   !if(allocated(PhosBot  ))
          deallocate(PAlk       )   !if(allocated(PAlk     ))
          deallocate(SiAlk      )   !if(allocated(SiAlk    ))
          deallocate(FREEtoTOT  )   !if(allocated(FREEtoTOT))
          deallocate(Hfree      )   !if(allocated(Hfree    ))
          deallocate(HSO4       )   !if(allocated(HSO4     ))
          deallocate(HF         )   !if(allocated(HF       ))
             
    end subroutine CalculateTAfromTCpH



    subroutine CalculatepHfromTCfCO2 &
               (TCi, fCO2i, pHctemp, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculatepHfromTCfCO2, version 02.02, 11-12-96, written by Ernie Lewis.
        ! ' Inputs: TC, fCO2, K0, K1, K2
        ! ' Output: pH
        ! ' This calculates pH from TC and fCO2 using K0, K1, and K2 by solving the
        ! '       quadratic in H: fCO2.*K0 = TC.*H.*H./(K1.*H + H.*H + K1.*K2).
        ! ' if there is not a real root, then pH is returned as missingn.
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCi
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: fCO2i
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHctemp

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: RR
        real(kind = DBL_PREC), allocatable, dimension(:) :: Discr
        real(kind = DBL_PREC), allocatable, dimension(:) :: H

        if(.not.(allocated(RR     ))) allocate(RR     (ntps))
        if(.not.(allocated(Discr  ))) allocate(Discr  (ntps))
        if(.not.(allocated(H      ))) allocate(H      (ntps))
        if(.not.(allocated(pHctemp))) allocate(pHctemp(ntps))
        
        RR      = K0 * fCO2i / TCi
        Discr   = (K1 * RR) * (K1 * RR) + 4.0D0 * (1 - RR) * (K1 * K2 *RR)
        H       = 0.5D0 * (K1 * RR + sqrt(Discr)) / (1 - RR)
        pHctemp = log(H) / log(0.1D0)
        
!    Deallocate auxilary arrays. Petras        
        deallocate(RR     )   !if(allocated(RR     ))
        deallocate(Discr  )   !if(allocated(Discr  ))
        deallocate(H      )   !if(allocated(H      ))
    end subroutine CalculatepHfromTCfCO2



    subroutine CalculateTCfrompHfCO2 &
               (pHi, fCO2i, TCctemp, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculateTCfrompHfCO2, version 01.02, 12-13-96, written by Ernie Lewis.
        ! ' Inputs: pH, fCO2, K0, K1, K2
        ! ' Output: TC
        ! ' This calculates TC from pH and fCO2, using K0, K1, and K2.
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pHi
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: fCO2i
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: TCctemp

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac        
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: H

        if(.not.(allocated(H      ))) allocate(H      (ntps))
        if(.not.(allocated(TCctemp))) allocate(TCctemp(ntps))
        
        H       = 10.0D0 ** (-pHi)
        TCctemp = K0 * fCO2i * (H * H + K1 * H + K1 * K2) / (H * H)
        
        !varargout{1}=TCctemp;
        
        deallocate(H)
    end subroutine CalculateTCfrompHfCO2


    subroutine RevelleFactor &
               (DEBUG, TAi, TCi, Revelle, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        integer, intent(in) :: DEBUG
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TAi
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCi
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: Revelle

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        ! global WhichKs;
        ! ' SUB RevelleFactor, version 01.03, 01-07-97, written by Ernie Lewis.
        ! ' Inputs: WhichKs%, TA, TC, K0, K(), T()
        ! ' Outputs: Revelle
        ! ' This calculates the Revelle factor (dfCO2/dTC)|TA/(fCO2/TC).
        ! ' It only makes sense to talk about it at pTot = 1 atm, but it is computed
        ! '       here at the given K(), which may be at pressure <> 1 atm. Care must
        ! '       thus be used to see if there is any validity to the number computed.

        real(kind = DBL_PREC), allocatable, dimension(:) :: TC0
        real(kind = DBL_PREC), allocatable, dimension(:) :: dTC
        real(kind = DBL_PREC), allocatable, dimension(:) :: pHc
        real(kind = DBL_PREC), allocatable, dimension(:) :: fCO2c
        real(kind = DBL_PREC), allocatable, dimension(:) :: fCO2plus
        real(kind = DBL_PREC), allocatable, dimension(:) :: fCO2minus
        real(kind = DBL_PREC), allocatable, dimension(:) :: TCi_INTERN

        if(.not.(allocated(TCi_INTERN))) allocate(TCi_INTERN(ntps))
        if(.not.(allocated(TC0)))        allocate(TC0       (ntps))
 
        TCi_INTERN = TCi
        TC0        = TCi_INTERN

        allocate(dTC(ntps), pHc(ntps), fCO2c(ntps))

        ! ' 1 umol/kg-SW
        dTC = 0.000001

        ! ' Find fCO2 at TA, TC + dTC
        TCi_INTERN = TC0 + dTC

        call CalculatepHfromTATC&
             (DEBUG, TAi, TCi_INTERN, pHc, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call CalculatefCO2fromTCpH &
             (TCi_INTERN, pHc, fCO2c, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)
        
        if(.not.(allocated(fCO2plus))) allocate(fCO2plus(ntps))
        fCO2plus = fCO2c

        ! ' Find fCO2 at TA, TC - dTC
        TCi_INTERN = TC0 - dTC

        call CalculatepHfromTATC &
             (DEBUG, TAi, TCi_INTERN, pHc, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        call CalculatefCO2fromTCpH &
             (TCi_INTERN, pHc, fCO2c, &
              ntps, &
              pHScale, WhichKs, WhoseKSO4, p1, p2, &
              TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
              RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
              fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
              FugFac, VPFac)

        if(.not.(allocated(fCO2minus))) allocate(fCO2minus(ntps))
        fCO2minus = fCO2c

        ! CalculateRevelleFactor:
        if(.not.(allocated(Revelle))) allocate(Revelle(ntps))        
        Revelle = (fCO2plus - fCO2minus) / dTC / ((fCO2plus + fCO2minus) / TCi_INTERN)
        !varargout{1}=Revelle;
        
!       Deallocate auxilary arrays. Petras        
                  deallocate(TC0       )   !if(allocated(TC0       ))
                  deallocate(dTC       )   !if(allocated(dTC       ))
                  deallocate(pHc       )   !if(allocated(pHc       ))
                  deallocate(fCO2c     )   !if(allocated(fCO2c     ))
                  deallocate(fCO2plus  )   !if(allocated(fCO2plus  ))
                  deallocate(fCO2minus )   !if(allocated(fCO2minus ))
                  deallocate(TCi_INTERN)   !if(allocated(TCi_INTERN))      
        
    end subroutine RevelleFactor


    subroutine CalculateAlkParts &
               (pHx, TCx, HCO3, CO3, BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB CalculateAlkParts, version 01.03, 10-10-97, written by Ernie Lewis.
        ! ' Inputs: pH, TC, K(), T()
        ! ' Outputs: HCO3, CO3, BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF
        ! ' This calculates the various contributions to the alkalinity.
        ! ' Though it is coded for H on the total pH scale, for the pH values occuring
        ! ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
        ! ' negligible) as long as the K Constants are on that scale.

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pHx
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TCx
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: HCO3
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: CO3
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: BAlk
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: OH
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: PAlk
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: SiAlk
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: Hfree
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: HSO4
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: HF

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosTop
        real(kind = DBL_PREC), allocatable, dimension(:) :: PhosBot
        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT

        
        if(.not.(allocated(H        ))) allocate(H        (ntps))
        if(.not.(allocated(HCO3     ))) allocate(HCO3     (ntps))
        if(.not.(allocated(CO3      ))) allocate(CO3      (ntps))
        if(.not.(allocated(BAlk     ))) allocate(BAlk     (ntps))
        if(.not.(allocated(OH       ))) allocate(OH       (ntps))
        if(.not.(allocated(PhosTop  ))) allocate(PhosTop  (ntps))
        if(.not.(allocated(PhosBot  ))) allocate(PhosBot  (ntps))
        if(.not.(allocated(PAlk     ))) allocate(PAlk     (ntps))
        if(.not.(allocated(SiAlk    ))) allocate(SiAlk    (ntps))
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
        if(.not.(allocated(Hfree    ))) allocate(Hfree    (ntps))
        if(.not.(allocated(HSO4     ))) allocate(HSO4     (ntps))
        if(.not.(allocated(HF       ))) allocate(HF       (ntps))
           
        H       = 10.0D0 ** (-pHx)
        HCO3    = TCx * K1 * H / (K1 * H + H * H + K1 * K2)
        CO3     = TCx * K1 * K2 / (K1 * H + H * H + K1 * K2)
        BAlk    = TB * KB / (KB + H)
        OH      = KW / H
        PhosTop = KP1 * KP2 * H + 2.0D0 * KP1 * KP2 * KP3 - H * H * H
        PhosBot =  H * H * H + KP1 * H * H + KP1 * KP2 * H + KP1 * KP2 * KP3
        PAlk    =  TP * PhosTop / PhosBot

        ! this is good to better than .0006*TP:
        ! PAlk = TP*(-H/(KP1+H) + KP2/(KP2+H) + KP3/(KP3+H))

        SiAlk     = TSi * KSi / (KSi + H)

        !' pH scale conversion factor
        FREEtoTOT = (1.0D0 + TS / KS)

        !' for H on the total scale
        Hfree     = H / FREEtoTOT

        !' since KS is on the free scale
        HSO4      = TS / (1.0D0 + KS / Hfree)

        !' since KF is on the free scale
        HF        = TF / (1.0D0 + KF / Hfree)

        !varargout{1} = HCO3; varargout{2} = CO3;   varargout{3} = BAlk;  varargout{4} = OH;
        !varargout{5} = PAlk; varargout{6} = SiAlk; varargout{7} = Hfree; varargout{8} = HSO4;
        !varargout{9} = HF;
        
!      Deallocate auxilary arrays. Petras        
               deallocate(H        )    !if(allocated(H        ))
               deallocate(PhosTop  )    !if(allocated(PhosTop  ))
               deallocate(PhosBot  )    !if(allocated(PhosBot  ))
               deallocate(FREEtoTOT)    !if(allocated(FREEtoTOT))
    end subroutine CalculateAlkParts


    subroutine CaSolubility &
               (Salt, TempC, Pdbar, TC, pH, OmegaCa, OmegaAr, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! '***********************************************************************
        ! ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
        ! ' Inputs: WhichKs%, Sal, TempCi, Pdbari, TCi, pHi, K1, K2
        ! ' Outputs: OmegaCa, OmegaAr
        ! ' This calculates omega, the solubility ratio, for calcite and aragonite.
        ! ' This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
        ! '       where Ksp is the solubility product (either KCa or KAr).
        ! '***********************************************************************
        ! ' These are from:
        ! ' Mucci, Alphonso, The solubility of calcite and aragonite in seawater
        ! '       at various salinities, temperatures, and one atmosphere total
        ! '       pressure, American Journal of Science 283:781-799, 1983.
        ! ' Ingle, S. E., Solubility of calcite in the ocean,
        ! '       Marine Chemistry 3:301-319, 1975,
        ! ' Millero, Frank, The thermodynamics of the carbonate system in seawater,
        ! '       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
        ! ' Ingle et al, The solubility of calcite in seawater at atmospheric pressure
        ! '       and 35%o salinity, Marine Chemistry 1:295-307, 1973.
        ! ' Berner, R. A., The solubility of calcite and aragonite in seawater in
        ! '       atmospheric pressure and 34.5%o salinity, American Journal of
        ! '       Science 276:713-730, 1976.
        ! ' Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
        ! ' Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
        ! '       boric acid, and the pHi of seawater, Limnology and Oceanography
        ! '       13:403-417, 1968.
        ! '***********************************************************************

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: Salt
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TempC
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: Pdbar
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: TC
        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pH
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: OmegaCa
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: OmegaAr

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: Ca
        real(kind = DBL_PREC), allocatable, dimension(:) :: Ar
        real(kind = DBL_PREC), allocatable, dimension(:) :: KCa
        real(kind = DBL_PREC), allocatable, dimension(:) :: KAr

        real(kind = DBL_PREC), allocatable, dimension(:) :: logKCa
        real(kind = DBL_PREC), allocatable, dimension(:) :: logKAr
        real(kind = DBL_PREC), allocatable, dimension(:) :: deltaVKCa
        real(kind = DBL_PREC), allocatable, dimension(:) :: KappaKCa
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKCafac
        real(kind = DBL_PREC), allocatable, dimension(:) :: deltaVKAr
        real(kind = DBL_PREC), allocatable, dimension(:) :: KappaKAr
        real(kind = DBL_PREC), allocatable, dimension(:) :: lnKArfac

        real(kind = DBL_PREC), allocatable, dimension(:) :: H
        real(kind = DBL_PREC), allocatable, dimension(:) :: CO3

        allocate(Ca(ntps)      , Ar(ntps)       , KCa(ntps )     , KAr(ntps)     , &
                 logKCa(ntps)  , logKAr(ntps)   , deltaVKCa(ntps), KappaKCa(ntps), &
                 lnKCafac(ntps), deltaVKAr(ntps), KappaKAr(ntps) , lnKArfac(ntps))

        where((WhichKs.ne.6).and.(WhichKs.ne.7))
            ! (below here, F isn't used, since almost always all rows match the above criterium,
            !  in all other cases the rows will be overwritten later on).

            ! CalculateCa:
            ! '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
            ! '       this is .010285.*Sali./35
            Ca = 0.02128D0 / 40.087D0 * (Salt / 1.80655D0)    ! ' in mol/kg-SW

            ! CalciteSolubility:
            ! '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
            logKCa = -171.9065D0 - 0.077993D0 * TempK + 2839.319D0 / TempK
            logKCa = logKCa + 71.595D0 * logTempK / log(10.0D0)
            logKCa = logKCa + (-0.77712D0 + 0.0028426D0 * TempK + 178.34D0 / TempK) * sqrSal
            logKCa = logKCa - 0.07711D0 * Salt + 0.0041249D0 * sqrSal * Salt

            ! '       sd fit = .01 (for Salt part, not part independent of Salt)
            KCa = 10.0D0 ** (logKCa) ! ' this is in (mol/kg-SW)^2

            ! AragoniteSolubility:
            ! '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
            logKAr = -171.945D0 - 0.077993D0 * TempK + 2903.293D0 /TempK
            logKAr = logKAr + 71.595D0 * logTempK / log(10.0D0)
            logKAr = logKAr + (-0.068393D0 + 0.0017276D0 * TempK + 88.135D0 /TempK) * sqrSal
            logKAr = logKAr - 0.10018D0 * Salt + 0.0059415D0 * sqrSal * Salt

            ! '       sd fit = .009 (for Salt part, not part independent of Salt)
            ! ' this is in (mol/kg-SW)^2
            KAr = 10.0D0 ** (logKAr)

            ! PressureCorrectionForCalcite:
            ! '       Ingle, Marine Chemistry 3:301-319, 1975
            ! '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
            ! '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
            deltaVKCa = -48.76D0  + 0.5304D0 * TempC
            KappaKCa  = (-11.76D0 + 0.3692D0 * TempC) / 1000.0D0
            lnKCafac  = (-deltaVKCa + 0.5D0 * KappaKCa * Pbar) * Pbar /RT
            KCa       = KCa * exp(lnKCafac)

            ! PressureCorrectionForAragonite:
            ! '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
            ! '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
            ! '       and 10^3 for Kappa factor)
            deltaVKAr = deltaVKCa + 2.8D0
            KappaKAr  = KappaKCa
            lnKArfac  = (-deltaVKAr + 0.5D0 * KappaKAr * Pbar) * Pbar / RT
            KAr       = KAr * exp(lnKArfac)
        end where

        where((WhichKs.eq.6).or.(WhichKs.eq.7))
            !
            ! *** CalculateCaforGEOSECS:
            ! Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
            ! (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
            Ca = 0.01026 * Salt/35
            ! Culkin gives Ca = (.0213./40.078).*(Sal./1.80655) in mol/kg-SW
            ! which corresponds to Ca = .01030.*Sal./35.

            ! *** CalculateKCaforGEOSECS:
            ! Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
            ! (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
            ! but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
            KCa = 0.0000001D0 * (-34.452D0 - 39.866 * Salt ** (1.0D0/3.0D0) + &
                 110.21D0 * log(Salt) / log(10.0D0) - 0.0000075752D0 * TempK * Tempk)
            ! this is in (mol/kg-SW)^2

            ! *** CalculateKArforGEOSECS:
            ! Berner, R. A., American Journal of Science 276:713-730, 1976:
            ! (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
            KAr = 1.45 * KCa ! this is in (mol/kg-SW)^2
            ! Berner (p. 722) states that he uses 1.48.
            ! It appears that 1.45 was used in the GEOSECS calculations

            ! *** CalculatePressureEffectsOnKCaKArGEOSECS:
            ! Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
            ! (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
            ! but their paper is not even on this topic).
            ! The fits appears to be new in the GEOSECS report.
            ! I can't find them anywhere else.
            KCa = KCa * exp((36D0   - 0.2D0  * TempC) * Pbar / RT)
            KAr = KAr * exp((33.3D0 - 0.22D0 * TempC) * Pbar / RT)
        end where

        ! CalculateOmegasHere:
        if(.not.(allocated(H      ))) allocate(H      (ntps))
        if(.not.(allocated(CO3    ))) allocate(CO3    (ntps))
        if(.not.(allocated(OmegaCa))) allocate(OmegaCa(ntps))
        if(.not.(allocated(OmegaAr))) allocate(OmegaAr(ntps))
        
        H       = 10.0D0 ** (-pH)
        CO3     = TC * K1 *K2 / (K1 * H + H * H + K1 * K2)
        OmegaCa = CO3 * Ca / KCa
        OmegaAr = CO3 * Ca / KAr
        !varargout{1} = CO3.*Ca./KCa; % OmegaCa, dimensionless
        !varargout{2} = CO3.*Ca./KAr; % OmegaAr, dimensionless
        
!        Deallocate auxilary arrays. Petras        
           deallocate(Ca       )    !if(allocated(Ca       ))
           deallocate(Ar       )    !if(allocated(Ar       ))
           deallocate(KCa      )    !if(allocated(KCa      ))
           deallocate(KAr      )    !if(allocated(KAr      ))
           deallocate(logKCa   )    !if(allocated(logKCa   ))
           deallocate(logKAr   )    !if(allocated(logKAr   ))
           deallocate(deltaVKCa)    !if(allocated(deltaVKCa))
           deallocate(KappaKCa )    !if(allocated(KappaKCa ))
           deallocate(lnKCafac )    !if(allocated(lnKCafac ))
           deallocate(deltaVKAr)    !if(allocated(deltaVKAr))
           deallocate(KappaKAr )    !if(allocated(KappaKAr ))
           deallocate(lnKArfac )    !if(allocated(lnKArfac ))
           deallocate(H        )    !if(allocated(H        ))
           deallocate(CO3      )    !if(allocated(CO3      ))  
    end subroutine CaSolubility


    subroutine FindpHOnAllScales &
               (pH, pHtot, pHsws, pHfree, pHNBS, &
                ntps, &
                pHScale, WhichKs, WhoseKSO4, p1, p2, &
                TempCi, TempCo, Pdbari, Pdbaro, Sal, sqrSal, TP, TSi, &
                RGasConstant, TB, TF, TS, TempK, RT, logTempK, Pbar, TempK100, lnK0, K0, KS, &
                fH, K1, K2, KW, KB, KF, KP1, KP2, KP3, KSi, &
                FugFac, VPFac)

        ! ' SUB FindpHOnAllScales, version 01.02, 01-08-97, written by Ernie Lewis.
        ! ' Inputs: pHScale%, pH, K(), T(), fH
        ! ' Outputs: pHNBS, pHfree, pHTot, pHSWS
        ! ' This takes the pH on the given scale and finds the pH on all scales.
        !  TS = T(3); TF = T(2);
        !  KS = K(6); KF = K(5);% 'these are at the given T, S, P

        real(kind = DBL_PREC), intent(in)   , allocatable, dimension(:) :: pH
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHtot
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHsws
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHfree
        real(kind = DBL_PREC), intent(inout), allocatable, dimension(:) :: pHNBS

        ! FORMER GLOBAL VARIABLES IN COSYS
        integer, intent(inout) :: ntps

        integer, intent(inout), dimension(ntps) :: pHScale
        integer, intent(inout), dimension(ntps) :: WhichKs
        integer, intent(inout), dimension(ntps) :: WhoseKSO4
        integer, intent(inout), dimension(ntps) :: p1
        integer, intent(inout), dimension(ntps) :: p2

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCi
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempCo
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbari
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pdbaro
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Sal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: sqrSal
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TP
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TSi

        real(kind = DBL_PREC) :: RGasConstant

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: RT
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: logTempK
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: Pbar
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: TempK100
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: lnK0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K0
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KS

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: fH
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: K2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KW
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KB
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KF
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP1
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP2
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KP3
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: KSi

        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: FugFac
        real(kind = DBL_PREC), intent(inout), dimension(ntps) :: VPFac
        ! END OF GLOBAL VARIABLES IN COSYS

        real(kind = DBL_PREC), allocatable, dimension(:) :: FREEtoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: SWStoTOT
        real(kind = DBL_PREC), allocatable, dimension(:) :: factor

        ! ' pH scale conversion factor
        if(.not.(allocated(FREEtoTOT))) allocate(FREEtoTOT(ntps))
        FREEtoTOT = (1.0D0 + TS / KS)

        ! ' pH scale conversion factor
        if(.not.(allocated(SWStoTOT))) allocate(SWStoTOT(ntps))
        SWStoTOT = (1.0D0 + TS / KS) / (1.0D0 + TS/KS + TF/KF)

        allocate(factor(ntps))

        !'"pHtot"
        where(pHScale==1)
            factor = 0.0D0
        end where

        ! '"pHsws"
        where(pHScale==2)
            factor = -log(SWStoTOT) / log(0.1D0)
        end where

        ! '"pHfree"
        where(pHScale==3)
            factor = -log(FREEtoTOT) / log(0.1D0)
        end where

        !'"pHNBS"
        where(pHScale==4)
            factor = -log(SWStoTOT) / log(0.1D0) + log(fH) / log(0.1D0)
        end where

        if(.not.(allocated(pHtot ))) allocate(pHtot (ntps))
        if(.not.(allocated(pHNBS ))) allocate(pHNBS (ntps))
        if(.not.(allocated(pHfree))) allocate(pHfree(ntps))
        if(.not.(allocated(pHsws ))) allocate(pHsws (ntps))
        
        pHtot = pH - factor ! ' pH comes into this sub on the given scale
        pHNBS = pHtot - log(SWStoTOT)  / log(0.1D0) + log(fH) / log(0.1D0)
        pHfree= pHtot - log(FREEtoTOT) / log(0.1D0)
        pHsws = pHtot - log(SWStoTOT)  / log(0.1D0)
        !varargout{1}=pHtot;
        !varargout{2}=pHsws;
        !varargout{3}=pHfree;
        !varargout{4}=pHNBS;
        
!       Deallocate auxilary arrays. Petras        
          deallocate(FREEtoTOT)   !if(allocated(FREEtoTOT))
          deallocate(SWStoTOT )   !if(allocated(SWStoTOT ))
          deallocate(factor   )   !if(allocated(factor   ))
        
    end subroutine FindpHOnAllScales

 end module CO2SYS_CDIAC

!******************************************************************
!******************************************************************
!******************************************************************
