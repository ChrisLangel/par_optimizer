      MODULE NORM_MOD
! 
      IMPLICIT NONE
      REAL(KIND=8) :: INTNORM,PTNORM 
!      
      CONTAINS
!         Function that computes norm of integral between a desired
!         value and the value of the current function with the current
!         parameter vector                         
          SUBROUTINE COMPUTEINTNORM
          USE DIMEN_MOD
          USE GRID_MOD
          USE CASE_MOD
!
          IMPLICIT NONE
          INTEGER :: I,J,K,N,G,P,IG,JG
          INTEGER :: ISHIFT,JSHIFT,NUMPAR,SIZEP,JIND,IND 
          REAL(KIND=8) :: INTAR,INTARHT,INTNORMC,AVG 
          INTEGER,DIMENSION(:),ALLOCATABLE :: JINDS
!
!
          INTNORM = 0.0 
! 
          DO N = 1,NLCASE
!         Need to do a bit of a coordinate transformation, the grid
!         file is indexed slightly different than the case
!        
          G      = GNUMS(N)          ! get the grid number 
          ISHIFT = ISA(N) - GSI(G)   ! look at grid starting index
          JSHIFT = JSA(N) - GSJ(G) 
!
!         ########################################################
!         Figure out the indicies of the domain partitions, thinking it
!         would be a good idea to pass info about this in cases file
          NUMPAR = 10                 ! Hard code for now
!         
          IF ( (JTA(N)-1) < NUMPAR ) THEN
            NUMPAR = JTA(N)-1 
          END IF
!
          IF (.NOT. ALLOCATED(JINDS) ) THEN
             ALLOCATE( JINDS((NUMPAR+2)) ) 
          END IF
          SIZEP    = (JTA(N)-1)/NUMPAR 
          JINDS(1) = 1  
          JINDS(NUMPAR+2) = JTA(N) - 1  
          DO K = 2,NUMPAR+1
             JINDS(K) = SIZEP*(K-1) 
          END DO 
!         ###########################################################       
!          
          INTNORMC = 0.0  ! The integral norm for this case                
          DO P = 1,NUMPAR+1  
!            Initialize sums 
             INTAR   = 0.0
             INTARHT = 0.0 
! 
             DO J = JINDS(P),JINDS(P+1)
                !WRITE(*,*) JINDS(P),JINDS(P+1)
                DO I = 1,ITA(N)-1 
                   IG      = I + ISHIFT
                   JG      = J + JSHIFT
                   AVG = (CASEVARS(N)%FAR(I,J)+CASEVARS(N)%FAR(I+1,J) &
     &           + CASEVARS(N)%FAR(I,J+1)+CASEVARS(N)%FAR(I+1,J+1))/4.0D0    
                   INTAR   = INTAR   + GRIDIN(G)%VOL(IG,JG)*AVG
!
                   AVG =(CASEVARS(N)%FHAT(I,J)+CASEVARS(N)%FHAT(I+1,J) &
     &           + CASEVARS(N)%FHAT(I,J+1)+CASEVARS(N)%FHAT(I+1,J+1))/4.0D0
                   INTARHT = INTARHT + GRIDIN(G)%VOL(IG,JG)*AVG
                END DO
             END DO
             INTNORMC = INTNORMC + ABS(INTARHT - INTAR) 
          END DO 
!          WRITE(*,*) 'The intnorm is',INTNORMC
          DEALLOCATE( JINDS ) 
!         
          INTNORM = INTNORM + INTNORMC
          END DO ! end case loop 
!
          END SUBROUTINE COMPUTEINTNORM

!        ---------------------------------------------------------------
!         Function that computes pointwise norm between a desired
!         value and the value of the current function with the current
!         parameter vector                         
!        ---------------------------------------------------------------
          SUBROUTINE COMPUTENORM
          USE DIMEN_MOD
!          USE GRID_MOD
          USE CASE_MOD
!
          IMPLICIT NONE
          INTEGER :: I,J,N
!
          PTNORM = 0.0
!
          DO N = 1,NLCASE
!
          DO J = 1,JTA(N)
             DO I = 1,ITA(N) 
                PTNORM = PTNORM + & 
     &          ABS(CASEVARS(N)%FAR(I,J) - CASEVARS(N)%FHAT(I,J))
!                WRITE(*,*) CASEVARS(N)%FAR(I,J),CASEVARS(N)%AR(I,J)
             END DO
          END DO 

          END DO 
!        
          END SUBROUTINE COMPUTENORM
!
      END MODULE NORM_MOD


