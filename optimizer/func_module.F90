      MODULE FUNC_MOD
! 
      IMPLICIT NONE
!      
      CONTAINS
!
!         This computes is the actual function we are trying to optimize,
!         routine will go through the case vars and produce the derived 
!         quantity F_Ar using known variables
!              
          SUBROUTINE COMPUTEFUNC 
          USE DIMEN_MOD
          USE PARAM_MOD 
          USE CASE_MOD
!
          IMPLICIT NONE
          INTEGER :: I,J,N
          REAL(KIND=8) :: FVAL 
!
          DO N = 1,NLCASE
             DO J = 1,JTA(N)
                DO I = 1,ITA(N)
                CALL FUNC(FVAL,CASEVARS(N)%LBTH(I,J), &
     &          CASEVARS(N)%AR(I,J),CASEVARS(N)%FTH(I,J) )
!                
                CASEVARS(N)%FAR(I,J) = FVAL  
                END DO
             END DO
          END DO 
!
          END SUBROUTINE COMPUTEFUNC

!        ---------------------------------------------------------------
!         ******* Function we are actually optimizing *************
!        ---------------------------------------------------------------
          SUBROUTINE FUNC(FVAL,LBT,AR,FTH)

!          USE DIMEN_MOD 
          USE PARAM_MOD
!
          IMPLICIT NONE
          REAL(KIND=8), INTENT(IN) :: LBT,AR,FTH 
          REAL(KIND=8), INTENT(OUT) :: FVAL  
          INTEGER :: I,J,N
!
!
          IF (LBT .GT. 1.0D0) THEN
             FVAL = (PMVEC(1)+PMVEC(5)*AR)*FTH*(AR**PMVEC(3)) + &
     &              PMVEC(2)*((1.0D0/LBT)**PMVEC(4))*(1-FTH)*AR 
          ELSE  
             FVAL = (PMVEC(1)+PMVEC(5)*AR)*FTH*(AR**PMVEC(3)) + &
     &              PMVEC(2)*(1.0D0/LBT)*(1-FTH)*AR
          END IF  
!        
          END SUBROUTINE FUNC
!
!        ---------------------------------------------------------------
!         ******* End of function we are actually optimizing *************
!        ---------------------------------------------------------------

!
          SUBROUTINE WRITEDIFF
          USE DIMEN_MOD
          IMPLICIT NONE
          INTEGER :: N,I,J,K 
          INTEGER :: IUNIT 
          CHARACTER(80) :: outname 
          REAL(KIND=8) :: ftemp 
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: DIFFAR
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: TVARS
!           
          DO N = 1,NLCASE
!            Allocate temporary array that will store info from file   
                     
             K = 4 ! Number of different vars in array #####   
             IF( .NOT. ALLOCATED(TVARS) ) THEN
                ALLOCATE( TVARS(GIT(N),GJT(N),K) ) 
             END IF 
!             
             IUNIT = (RK*10)+NPROCS*10 
!             WRITE(*,*) 'Look at casename', casename,IUNIT 
             OPEN(FILE=lcnames(N),UNIT=IUNIT,STATUS='OLD',&
     &             ACTION='READ',FORM='UNFORMATTED')
             REWIND(IUNIT) 
             READ(IUNIT)  ! Read lines from header 
             READ(IUNIT)
             READ(IUNIT)
!            Read actual array that contains vars             
             READ(IUNIT) TVARS 
!
             CLOSE(IUNIT) 
             IF ( .NOT. ALLOCATED(DIFFAR) ) THEN 
                ALLOCATE( DIFFAR(GIT(N),GJT(N)) )
             END IF  
!            Go through and compute the difference at all points
             DO J = 1, GJT(N) 
                DO I = 1,GIT(N) 
                   CALL FUNC(ftemp,TVARS(I,J,1),TVARS(I,J,2), &
     &                       TVARS(I,J,3)  ) 
                   
                   DIFFAR(I,J) = ftemp - TVARS(I,J,4) 
                END DO 
             END DO 
             outname = trim(lcnames(N)) // "_DIFF"

!            Write out the difference array to file
             OPEN(IUNIT,FILE=outname,STATUS='REPLACE', &
     &           FORM='UNFORMATTED')
             WRITE(IUNIT) GIT(N),GJT(N)
             WRITE(IUNIT) DIFFAR 
!             
             CLOSE(IUNIT) 
             DEALLOCATE( DIFFAR ) 
          END DO
          END SUBROUTINE WRITEDIFF 
! 
!
!
      END MODULE FUNC_MOD


