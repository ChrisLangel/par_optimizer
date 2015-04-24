      MODULE VOL_MOD
!
      IMPLICIT NONE
!      
!  
      CONTAINS 
         SUBROUTINE COMPUTEVOL 
         USE GRID_MOD
         USE DIMEN_MOD  
         IMPLICIT NONE
         INTEGER :: I,J,K,N,MAXI,MAXJ,ITEMP,JTEMP
         INTEGER :: IS,IE,JS,JE
         REAL(KIND=8) :: HALF
!
!
         HALF = 0.05D0 
         DO K = 1,NLGRIDS
            IE = GEI(K) - GSI(K) 
            JE = GEJ(K) - GSJ(K)       
            DO J = 1,JE
               DO I = 1,IE
                  GRIDIN(K)%VOL(I,J) = HALF* & 
              ( ( GRIDIN(K)%X(I,J)*    GRIDIN(K)%Y(I+1,J)   -   &
                  GRIDIN(K)%X(I+1,J)*  GRIDIN(K)%Y(I,J)     ) + &     
                ( GRIDIN(K)%X(I+1,J)*  GRIDIN(K)%Y(I+1,J+1) -   &
                  GRIDIN(K)%X(I+1,J+1)*GRIDIN(K)%Y(I+1,J)   ) + &    
                ( GRIDIN(K)%X(I+1,J+1)*GRIDIN(K)%Y(I,J+1)   -   &
                  GRIDIN(K)%X(I,J+1)*  GRIDIN(K)%Y(I+1,J+1) ) + &    
                ( GRIDIN(K)%X(I,J+1)*  GRIDIN(K)%Y(I,J)     -   &
                  GRIDIN(K)%X(I,J)*    GRIDIN(K)%Y(I,J+1)   )   )  
               END DO 
            END DO 
         END DO   
         END SUBROUTINE
!
      END MODULE VOL_MOD


