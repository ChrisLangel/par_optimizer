      MODULE OPT_MOD
! 
      USE DIMEN_MOD
      USE PARAM_MOD      
      IMPLICIT NONE
      INTEGER :: MAXITS
      REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: GRADVEC,STEPDIR,&
     &                                         PMVECK,GRADVECK  
      REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Hk
      REAL(KIND=8) :: DIFFSTEP,STEPFRAC,INTNCN,PTNCON,TOL 
      CHARACTER(16) :: otype                
!      
      CONTAINS
!              
!        ---------------------------------------------------------------
!         This computes the objective function using the norm routines
!         in NORM_MOD
!        --------------------------------------------------------------- 
          SUBROUTINE OBJFUNC 
          USE PARAM_MOD
          USE FUNC_MOD
          USE NORM_MOD 
!
          IMPLICIT NONE

!          INTNCN = 1.0D0
!          PTNCON = 0.0010D0
!
          CALL COMPUTEFUNC 
          CALL COMPUTEINTNORM
          CALL COMPUTENORM 
!            
          TFVAL = INTNCN*INTNORM +PTNCON*PTNORM 
!
          END SUBROUTINE OBJFUNC
!        ---------------------------------------------------------------
!   

!        ---------------------------------------------------------------
!         This computes finite difference gradients with respect to the
!         vector of constants 
!        ---------------------------------------------------------------
!              
          SUBROUTINE GRADFUNC 
          USE DIMEN_MOD
          USE PARAM_MOD
!
          IMPLICIT NONE
          INTEGER :: N
          REAL(KIND=8) :: DF,DX,MAG 
          REAL(KIND=8), DIMENSION(NCON) :: TPMVEC 
          
!          DIFFSTEP = 0.00000000000001D0  
!         Establish initial function value, and store the initial
!         parameter vector
          CALL OBJFUNC
          LFVAL  = TFVAL 
          TPMVEC = PMVEC 
!          

!         Go through and perturb the each of the constants in the vector
!         of coefficients
!          
          DO N = 1,NCON
             PMVEC      = TPMVEC
             PMVEC(N)   = PMVEC(N) + DIFFSTEP*PMVEC(N)               
             CALL OBJFUNC
!             WRITE(*,*) 'Tfval is: ', TFVAL,LFVAL
             DF         = TFVAL - LFVAL
             DX         = PMVEC(N) - TPMVEC(N)   
             GRADVEC(N) = DF/DX
          END DO              
!        
!         Normalize the gradient vector 
          MAG = 0.0
          DO N = 1,NCON
             MAG = MAG + GRADVEC(N)**2
          END DO 
          MAG = SQRT(MAG) 
          DO N = 1,NCON
             GRADVEC(N) = GRADVEC(N)/MAG
          END DO          
!
!         This will be overwritten if using a different method
          STEPDIR = GRADVEC     
!          
!         Make sure to reset parameter vector and objective function          
          PMVEC = TPMVEC
          CALL OBJFUNC

!          WRITE(*,*) 'Df and Dx',DF,DX,GRADVEC(1)  
!
          END SUBROUTINE GRADFUNC

!        ---------------------------------------------------------------
!        ---------------------------------------------------------------
!        This computes a line search based on the gradient vector until
!        the objective function reduced in value (1 proc version)  
!        ---------------------------------------------------------------
!              
          SUBROUTINE LINESEARCH
          USE DIMEN_MOD
          USE PARAM_MOD
          USE FUNC_MOD
          USE NORM_MOD 
!
          IMPLICIT NONE
          INTEGER :: N
          REAL(KIND=8) :: STEPLEN 
          REAL(KIND=8), DIMENSION(NCON) :: TPMVEC 
!
!         Record starting function value
          LFVAL  = TFVAL
          TPMVEC = PMVEC 
!         Try with a step length of one 
          STEPLEN = 1.0D0 
          DO N = 1,NCON
             PMVEC(N)   = TPMVEC(N) - STEPLEN*STEPDIR(N)
          ENDDO            
          CALL OBJFUNC
!         Backtracking line search 
          DO WHILE ( TFVAL > LFVAL ) 
          STEPLEN = STEPLEN*STEPFRAC
          DO N = 1,NCON
             PMVEC(N)   = TPMVEC(N) - STEPLEN*STEPDIR(N)
!             WRITE(*,*) 'The steplen is',STEPLEN 
          ENDDO            
          CALL OBJFUNC
!          WRITE(*,*) 'From line search', TFVAL,STEPLEN,GRADVEC(N)  
          END DO   

!          WRITE(*,*) 'The final value of Param, FVAL:', PMVEC(1),TFVAL 
          END SUBROUTINE LINESEARCH
!
!
!        ---------------------------------------------------------------
!        This computes a line search based on the gradient vector until
!        the objective function reduced in value (Multi proc version)  
!        ---------------------------------------------------------------
!              
          SUBROUTINE LINESEARCHP 
          USE DIMEN_MOD
          USE PARAM_MOD 
          USE FUNC_MOD
          USE NORM_MOD 
!
          IMPLICIT NONE
          INTEGER :: N,ITER,MAXITER
          REAL(KIND=8) :: STEPLEN 
          REAL(KIND=8), DIMENSION(NCON) :: TPMVEC 
!
!          MAXITER = 100
!         Record starting function value
          LFVAL  = TFVAL
          TPMVEC = PMVEC 
!         Try with a step length of one 
          STEPLEN = 1.0D0 
          DO N = 1,NCON
             PMVEC(N)   = TPMVEC(N) - STEPLEN*STEPDIR(N)
          ENDDO            
          CALL SUMOBJFUNC  
          ITER = 1
!         Backtracking line search 
          DO WHILE ( TFVAL > LFVAL ) 
          STEPLEN = STEPLEN*STEPFRAC
          DO N = 1,NCON
             PMVEC(N)   = TPMVEC(N) - STEPLEN*STEPDIR(N)
!             WRITE(*,*) 'Grad val,ITer,RK',GRADVEC(N),ITER,RK 
          ENDDO            
          CALL SUMOBJFUNC
          ITER = ITER + 1
!          WRITE(*,*) 'From line search', TFVAL,STEPLEN,GRADVEC(N)  
          END DO   

!          WRITE(*,*) 'The final value of Param, FVAL:',PMVEC(1),&
!     &                                       PMVEC(2),TFVAL 
          END SUBROUTINE LINESEARCHP 
!

!        ---------------------------------------------------------------
!         This gathers all the gradient information from the procs and
!         computes a global gradient vec
!        ---------------------------------------------------------------
!              
          SUBROUTINE SUMGRAD 
          USE DIMEN_MOD
          USE PARAM_MOD
!
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: I,N,TAGB,TAG,IERR
          INTEGER :: STATUS(MPI_STATUS_SIZE)
          REAL(KIND=8) :: MAG 
          REAL(KIND=8),DIMENSION(NPROCS,NCON) :: RECA     
!
          CALL GRADFUNC
!         
          IF (RK .ne. 0) THEN
            TAG  = RK*22
            TAGB = RK*37
            CALL MPI_SEND(GRADVEC,NCON,MPI_DOUBLE_PRECISION,&
     &                        0,TAG,MPI_COMM_WORLD,IERR)  
             
            CALL MPI_RECV(GRADVEC,NCON,MPI_DOUBLE_PRECISION,&
     &                        0,TAGB,MPI_COMM_WORLD,STATUS,IERR) 

!             WRITE(*,*) 'The summed Gradient vector:', GRADVEC,RK     
          ELSE 
             DO N = 1,NPROCS-1 
                TAG = N*22 
                CALL MPI_RECV(RECA(N,1:NCON),NCON,MPI_DOUBLE_PRECISION,&
     &                        N,TAG,MPI_COMM_WORLD,STATUS,IERR) 
             END DO
             RECA(NPROCS,1:NCON) = GRADVEC
             DO I =1,NCON
                DO N = 1,NPROCS
                   GRADVEC(I) = GRADVEC(I) + RECA(N,I) 
                END DO
             END DO

!            Normalize the gradient vector 
             MAG = 0.0
             DO I = 1,NCON
                MAG = MAG + GRADVEC(I)**2
             END DO 
             MAG = SQRT(MAG) 
             DO I = 1,NCON
                GRADVEC(I) = GRADVEC(I)/MAG
             END DO
!
!            Send to all other procs
             DO N = 1,NPROCS-1
                TAGB = N*37
                CALL MPI_SEND(GRADVEC,NCON,MPI_DOUBLE_PRECISION,&
     &                        N,TAGB,MPI_COMM_WORLD,IERR)  
             END DO 
          END IF 
          CALL SUMOBJFUNC 
!
!         Again this will be overriden if using alternate method          
          STEPDIR = GRADVEC 
!            
!
          END SUBROUTINE SUMGRAD
!         
!        ---------------------------------------------------------------
!         This gathers all procs contributions to the objective function 
!         and computes a global function val 
!        ---------------------------------------------------------------
          SUBROUTINE SUMOBJFUNC 
          USE DIMEN_MOD
          USE PARAM_MOD 
!         USE FUNC_MOD
!          USE NORM_MOD 
!
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: I,N,TAGB,TAG,IERR
          INTEGER :: STATUS(MPI_STATUS_SIZE)
          REAL(KIND=8),DIMENSION(NPROCS) :: FARY
!
          CALL OBJFUNC
!         
          IF (RK .ne. 0) THEN
            TAG  = RK*23
            TAGB = RK*31
            CALL MPI_SEND(TFVAL,1,MPI_DOUBLE_PRECISION,&
     &                    0,TAG,MPI_COMM_WORLD,IERR)  
!            WRITE(*,*) 'Sending this val:',TFVAL 
            CALL MPI_RECV(TFVAL,1,MPI_DOUBLE_PRECISION,&
     &                    0,TAGB,MPI_COMM_WORLD,STATUS,IERR) 

!             WRITE(*,*) 'The summed function value vector:',TFVAL,RK   
          ELSE
!            Recieve contribution from all other procs               
             DO N = 1,NPROCS-1 
                TAG = N*23
                CALL MPI_RECV(FARY(N),1,MPI_DOUBLE_PRECISION,&
     &                        N,TAG,MPI_COMM_WORLD,STATUS,IERR) 
             END DO
!            Insert root procs contribution into array            
             FARY(NPROCS) = TFVAL 
             TFVAL        = 0.0
             DO N = 1,NPROCS
                TFVAL = TFVAL + FARY(N)  
             END DO
!
!            Send to all other procs
             DO N = 1,NPROCS-1
                TAGB = N*31
                CALL MPI_SEND(TFVAL,1,MPI_DOUBLE_PRECISION,&
     &                        N,TAGB,MPI_COMM_WORLD,IERR)  
             END DO 
          END IF 
 
!            
!
          END SUBROUTINE SUMOBJFUNC


!        ---------------------------------------------------------------
!         This computes the search direction for a Quasi-Newton method
!         (Multi-Proc) 
!        ---------------------------------------------------------------
!              
          SUBROUTINE QUASINEWP
          USE DIMEN_MOD
          USE PARAM_MOD 
!
          IMPLICIT NONE
          INTEGER :: I,J,N
          REAL(KIND=8) :: RHOK,MAG
          REAL(KIND=8), DIMENSION(NCON) :: TPMVEC,PK,SK,YK
          REAL(KIND=8), DIMENSION(NCON,NCON) :: HnOT,HpONE,TMPL,TMPR, &
     &                                          TMPA  
!
          TMPL = 0.0D0
          TMPR = 0.0D0
          HnOT = 0.0D0 
!         Initialize H_o to the identity matrix
          DO I = 1,NCON
             TMPL(I,I) = 1.0D0
             TMPR(I,I) = 1.0D0
             HnOT(I,I) = 1.0D0
          END DO 
!          
          CALL SUMGRAD
!         Now we have the gradient vector with contributions from all
!         procs, as well as global function value stored locally
!
          LFVAL = TFVAL   ! Store the function value at the start  
          
          IF (OPTITER == 1) THEN 
             CALL DOTMV(HnOT,GRADVEC,NCON,PK)
             !PK      = MATMUL(HnOT,GRADVEC) 
             Hk      = HnOT             
          ELSE 
             DO I = 1,NCON
                SK(I) = PMVEC(I)   - PMVECK(I)
                YK(I) = GRADVEC(I) - GRADVECK(I)
             END DO   
             CALL DOTVV(SK,YK,NCON,RHOK)  
             IF (RHOK > 0) THEN
                RHOK  = 1.0D0/RHOK
             END IF
             DO I = 1,NCON
                DO J = 1,NCON
                   TMPL(I,J) = TMPL(I,J) - RHOK*SK(I)*YK(J)
                   TMPR(I,J) = TMPR(I,J) - RHOK*SK(J)*YK(I)
                   TMPA(I,J) = RHOK*SK(I)*SK(J)  
                END DO
             END DO
!             WRITE(*,*) 'From QN routine TMPL', TMPL(1,1),TMPL(1,2)
!             WRITE(*,*) '--------------------', TMPL(2,1),TMPL(2,2)
!             WRITE(*,*) 'From QN routine TMPR', TMPR(1,1),TMPR(1,2)
!             WRITE(*,*) '--------------------', TMPR(2,1),TMPR(2,2)
             HpONE   = MATMUL(MATMUL(TMPL,Hk),TMPR) 
             DO I = 1,NCON
                DO J = 1,NCON
                   HpONE(I,J) = HpONE(I,J) + TMPA(I,J) 
                END DO
             END DO  
             PK      = MATMUL(HpONE,GRADVEC) 
             Hk      = HpONE
          END IF
          STEPDIR  = PK
!         Normalize the step direction vector 
          MAG = 0.0
          DO N = 1,NCON
             MAG = MAG + STEPDIR(N)**2
          END DO 
          MAG = SQRT(MAG) 
          DO N = 1,NCON
             STEPDIR(N) = STEPDIR(N)/MAG
          END DO 
!          

          GRADVECK = GRADVEC
          PMVECK   = PMVEC          
!          IF (RK == 0) THEN 
!          WRITE(*,*) 'From QN routine RHOK', RHOK
!          WRITE(*,*) 'From QN routine PK', STEPDIR(1),STEPDIR(2)
!          WRITE(*,*) 'From QN routine GRAD', GRADVEC(1),GRADVEC(2)
!          
!          WRITE(*,*) 'From QN routine SK', SK(1),SK(2)
!          WRITE(*,*) 'From QN routine YK', YK(1),YK(2)
!          WRITE(*,*) 'From QN routine Hk', Hk(1,1),Hk(1,2)
!          WRITE(*,*) '------------------', Hk(2,1),Hk(2,2)
!          END IF 
!        
!         Make sure to reset parameter vector and objective function          
          !PMVEC = TPMVEC

          END SUBROUTINE QUASINEWP


!        ---------------------------------------------------------------
!         This computes the search direction for a Quasi-Newton method
!         (Single-Proc) 
!        ---------------------------------------------------------------
!              
          SUBROUTINE QUASINEW
          USE DIMEN_MOD
          USE PARAM_MOD 
!
          IMPLICIT NONE
          INTEGER :: I,J,N
          REAL(KIND=8) :: RHOK,MAG
          REAL(KIND=8), DIMENSION(NCON) :: TPMVEC,PK,SK,YK
          REAL(KIND=8), DIMENSION(NCON,NCON) :: HnOT,HpONE,TMPL,TMPR, &
     &                                          TMPA  
!
          TMPL = 0.0D0
          TMPR = 0.0D0
          HnOT = 0.0D0 
!         Initialize H_o to the identity matrix
          DO I = 1,NCON
             TMPL(I,I) = 1.0D0
             TMPR(I,I) = 1.0D0
             HnOT(I,I) = 1.0D0
          END DO 
!          
          CALL GRADFUNC 
!         Now we have the gradient vector with contributions from all
!         procs, as well as global function value stored locally
!
          LFVAL = TFVAL   ! Store the function value at the start  
          
          IF (OPTITER == 1) THEN 
             CALL DOTMV(HnOT,GRADVEC,NCON,PK)
             !PK      = MATMUL(HnOT,GRADVEC) 
             Hk      = HnOT             
          ELSE 
             DO I = 1,NCON
                SK(I) = PMVEC(I)   - PMVECK(I)
                YK(I) = GRADVEC(I) - GRADVECK(I)
             END DO   
             CALL DOTVV(SK,YK,NCON,RHOK)  
             IF (RHOK > 0) THEN
                RHOK  = 1.0D0/RHOK
             END IF
             DO I = 1,NCON
                DO J = 1,NCON
                   TMPL(I,J) = TMPL(I,J) - RHOK*SK(I)*YK(J)
                   TMPR(I,J) = TMPR(I,J) - RHOK*SK(J)*YK(I)
                   TMPA(I,J) = RHOK*SK(I)*SK(J)  
                END DO
             END DO
!             WRITE(*,*) 'From QN routine TMPL', TMPL(1,1),TMPL(1,2)
!             WRITE(*,*) '--------------------', TMPL(2,1),TMPL(2,2)
!             WRITE(*,*) 'From QN routine TMPR', TMPR(1,1),TMPR(1,2)
!             WRITE(*,*) '--------------------', TMPR(2,1),TMPR(2,2)
             HpONE   = MATMUL(MATMUL(TMPL,Hk),TMPR) 
             DO I = 1,NCON
                DO J = 1,NCON
                   HpONE(I,J) = HpONE(I,J) + TMPA(I,J) 
                END DO
             END DO  
             PK      = MATMUL(HpONE,GRADVEC) 
             Hk      = HpONE
          END IF
          STEPDIR  = PK
!         Normalize the step direction vector 
          MAG = 0.0
          DO N = 1,NCON
             MAG = MAG + STEPDIR(N)**2
          END DO 
          MAG = SQRT(MAG) 
          DO N = 1,NCON
             STEPDIR(N) = STEPDIR(N)/MAG
          END DO 
          GRADVECK = GRADVEC
          PMVECK   = PMVEC        
!          
          END SUBROUTINE QUASINEW 



!
!      Subroutines that compute vector dot product and matrix vector 
!       multiplication         
          SUBROUTINE DOTVV(X,Y,LENGTH,DOTP)
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: LENGTH
          REAL(KIND=8), DIMENSION(LENGTH), INTENT(IN) :: X,Y
          REAL(KIND=8), INTENT(OUT) :: DOTP
          INTEGER :: I
          DOTP = 0.0
          DO I = 1,LENGTH
             DOTP = DOTP + X(I)*Y(I)
          END DO
          END SUBROUTINE DOTVV
!
!
          SUBROUTINE DOTMV(X,Y,LENGTH,DOTP)
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: LENGTH
          REAL(KIND=8), DIMENSION(LENGTH), INTENT(IN) :: Y
          REAL(KIND=8), DIMENSION(LENGTH,LENGTH), INTENT(IN) :: X
          REAL(KIND=8), DIMENSION(LENGTH), INTENT(OUT) :: DOTP
          INTEGER :: I,J 
          DOTP = 0.0D0
          DO I = 1,LENGTH
             DO J = 1,LENGTH
                DOTP(I) = DOTP(I) + X(I,J)*Y(I)
             END DO 
          END DO 
          END SUBROUTINE DOTMV

      END MODULE OPT_MOD


