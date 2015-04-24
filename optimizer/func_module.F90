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
      END MODULE FUNC_MOD


