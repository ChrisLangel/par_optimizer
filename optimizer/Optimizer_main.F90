      PROGRAM OPTIMIZE
!
! 
      USE DIMEN_MOD
      USE PARAM_MOD  
      USE GRID_MOD
      USE VOL_MOD
      USE CASE_MOD
      USE FUNC_MOD 
      USE NORM_MOD
      USE OPT_MOD
!      USE SECPROJ_MOD
!      USE SECDERIV_MOD
!      USE TIMEMARCH_MOD
!      USE RESID_MOD
!      USE CONNECT_MOD
!      USE COM_MOD
!
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER :: I,J,G,K,IERR,RANK,RC,NUMPROC, &
     &            clock_max,clock_rate     
      REAL(KIND=8) :: wall_time,MRESD,clock_start,clock_end
      CHARACTER(80) :: residname,qname,cname
!
!
!    ########## MPI INIT (Including Domain Partitioning) ################
!     
      CALL MPI_INIT(IERR)
!     
      IF (IERR .NE. MPI_SUCCESS) THEN
         print *,'Error starting MPI program. Terminating.'
         CALL MPI_ABORT(MPI_COMM_WORLD, RC, IERR)
      END IF 
!
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
!
      RK = RANK
      NPROCS = NUMPROC
!  
!     Have root processor partition domain and write out conectivity files
      IF (RANK == 0) THEN 
         CALL DIVIDECASES 
      END IF 
! 
!    ####################################################################
!   
!     Initialization routines 
      CALL READ_INPUT  ! Read input to determine info about run
      WRITE(*,*) 'The paramvec: ',PMVEC(1),PMVEC(2),PMVEC(3) 
      CALL GETGRIDCORDS
      CALL COMPUTEVOL  ! Compute and store cell volumes
      CALL READCASES   ! Store variables from each case
!     
!
      clock_start = MPI_WTIME()
      OPTITER = 1
!      
!     Don't even bother with the MPI_COM stuff if only run on one
!     processor 
      IF (NPROCS == 1) THEN
      CALL GRADFUNC  
      WRITE(*,*) 'The starting value/ Gradvec:',TFVAL,GRADVEC(1),&
     &                                       GRADVEC(2),GRADVEC(3)
      IF (otype (1:2) .eq. 'QN') THEN
         DO WHILE (OPTITER < MAXITS .and. TFVAL > TOL) 
         CALL QUASINEW
         CALL LINESEARCH
         OPTITER = OPTITER + 1
         WRITE(*,10) OPTITER,PMVEC(1),PMVEC(2),TFVAL
         END DO 
      ELSE  
         DO WHILE (OPTITER < MAXITS .and. TFVAL > TOL) 
         CALL GRADFUNC
         CALL LINESEARCH
         OPTITER = OPTITER + 1
         WRITE(*,10) OPTITER,PMVEC(1),PMVEC(2),TFVAL
         END DO 
      END IF 
      ELSE !################## MULTI-PROC #############################
      CALL SUMGRAD 
      WRITE(*,*) 'The starting value/ Gradvec:',TFVAL,GRADVEC(1),&
     &                                     GRADVEC(2),GRADVEC(3) 
      IF (otype(1:2) .eq. 'QN') THEN       
         DO WHILE (OPTITER < MAXITS .and. TFVAL > TOL)  
         CALL QUASINEWP  
         CALL LINESEARCHP
         OPTITER = OPTITER + 1
         DO I = 1,1-RK
         WRITE(*,10) OPTITER,PMVEC(1),PMVEC(2),TFVAL
         END DO 
         END DO
      ELSE  
         DO WHILE (OPTITER < MAXITS .and. TFVAL > TOL)  
         CALL SUMGRAD
         CALL LINESEARCHP
         OPTITER = OPTITER + 1
         DO I = 1,1-RK
         WRITE(*,10) OPTITER,PMVEC(1),PMVEC(2),TFVAL
         END DO 
         END DO 
      END IF  
      END IF 

 10   FORMAT('The value of Params at this Iter: ',I3,2F8.4,F16.12)

!    ---------------------------------------------------------------
!
      WRITE(*,*) 'The value of Param, FVAL at this iter:',I,PMVEC(1),&
     &             PMVEC(2),TFVAL
!      
      WRITE(*,*) 'The final value of Param, FVAL:',OPTITER,PMVEC(1),&
     &             PMVEC(2),PMVEC(3),PMVEC(4),TFVAL
!
!              
      clock_end =  MPI_WTIME()
      wall_time= clock_end-clock_start
      IF (RANK == 0) THEN 
         PRINT*,'Solver wall clocktime (seconds)',wall_time
      END IF 
        
!              
      IF (RK == 0) THEN 
         CLOSE (UNIT=999, STATUS="DELETE")
      END IF

      CALL MPI_FINALIZE(IERR)
      
!
      STOP
      END
!
