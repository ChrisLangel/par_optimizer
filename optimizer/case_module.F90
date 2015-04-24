      MODULE CASE_MOD
! 
      IMPLICIT NONE
! 
      TYPE VARS
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: LBTH,AR,FTH,FHAT,FAR
      END TYPE VARS
!
      TYPE(VARS), DIMENSION(:), ALLOCATABLE, SAVE :: CASEVARS  
!
!      
      CONTAINS
!         Function that reads the p3d var files and inserts into CASEVAR   
!              
          SUBROUTINE READCASES
          USE DIMEN_MOD
!
          IMPLICIT NONE
          INTEGER :: I,J,K,N,M,G,NI,NK,ITEMP,JTEMP,ICT,JCT,MAXI,MAXJ
          INTEGER :: IUNIT 
          CHARACTER(80) :: casename
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: TVARS
!
          IF (.NOT. ALLOCATED(CASEVARS) ) THEN
             ALLOCATE( CASEVARS(NLCASE) )
          END IF
!
!         Case arrays will only be the size of the region that
!         is being looked at, but need to make sure max is allocated 
!        
          MAXI = 0
          MAXJ = 0
          DO N = 1,NLCASE
             ITEMP = IEA(N)-ISA(N)+2 
             JTEMP = JEA(N)-JSA(N)+2
             IF (ITEMP > MAXI) THEN
                MAXI = ITEMP
             END IF
             IF (JTEMP > MAXJ) THEN
                MAXJ = JTEMP
             END IF
          END DO 
!          
          ALLOCATE( ITA(NLCASE),JTA(NLCASE) )           
!
!         Look through all the cases on this proc
!
          DO N = 1,NLCASE
!            Allocate temporary array that will store info from file   
                     
             K = 4 ! Number of different vars in array #####   
             IF( .NOT. ALLOCATED(TVARS) ) THEN
                ALLOCATE( TVARS(GIT(N),GJT(N),K) ) 
             END IF 
!
!            Allocate arrays to max size on this proc
             ALLOCATE( CASEVARS(N)%LBTH(MAXI,MAXJ) )          
             ALLOCATE( CASEVARS(N)%AR(MAXI,MAXJ)   )          
             ALLOCATE( CASEVARS(N)%FTH(MAXI,MAXJ)  )          
             ALLOCATE( CASEVARS(N)%FHAT(MAXI,MAXJ)  )          
             ALLOCATE( CASEVARS(N)%FAR(MAXI,MAXJ)  )          
!             
             casename = lcnames(N)
             IUNIT = (RK*10)+NPROCS*10 
!             WRITE(*,*) 'Look at casename', casename,IUNIT 
             OPEN(FILE=casename,UNIT=IUNIT,STATUS='OLD',&
     &            ACTION='READ',FORM='UNFORMATTED')
             REWIND(IUNIT) 
             READ(IUNIT)  ! Read lines from header 
             READ(IUNIT)
             READ(IUNIT)
!            Read actual array that contains vars             
             READ(IUNIT) TVARS 
!            
             JCT = 1
             DO J = JSA(N),JEA(N)
                ICT = 1
                DO I = ISA(N),IEA(N)
                   CASEVARS(N)%LBTH(ICT,JCT) = TVARS(I,J,1)
                   CASEVARS(N)%AR(ICT,JCT)   = TVARS(I,J,3)
                   CASEVARS(N)%FTH(ICT,JCT)  = TVARS(I,J,2)
                   CASEVARS(N)%FHAT(ICT,JCT)  = TVARS(I,J,4)
!                   WRITE(*,*) CASEVARS(N)%LBTH(ICT,JCT), &
!     &             CASEVARS(N)%AR(ICT,JCT),CASEVARS(N)%FTH(ICT,JCT)
                   ICT = ICT + 1
                END DO
                JCT = JCT + 1
             END DO
!
!             
             CLOSE(IUNIT)
             DEALLOCATE(TVARS)
             ITA(N) = IEA(N) - ISA(N) + 1
             JTA(N) = JEA(N) - JSA(N) + 1 
          END DO ! end grid loop 
!

          END SUBROUTINE READCASES 

!           
! 
!
      END MODULE CASE_MOD


