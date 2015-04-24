      MODULE GRID_MOD
! 
      IMPLICIT NONE
! 
      TYPE GRID
      REAL(KIND=8),  DIMENSION(:,:), ALLOCATABLE :: X,Y,VOL
      END TYPE GRID
!
      TYPE(GRID), DIMENSION(:), ALLOCATABLE, SAVE :: GRIDIN 
      CHARACTER(80), DIMENSION(:), ALLOCATABLE,SAVE :: lgnames ! watch 80
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: MGRID,NGRID,NIA,NJA
      INTEGER, SAVE :: NLGRIDS,SMGRID,SNGRID
!
!      
      CONTAINS
!         Function that reads the p3d grid files and inserts into GRIDIN              
          SUBROUTINE READGRIDS
          USE DIMEN_MOD
!
          IMPLICIT NONE
          INTEGER :: I,J,K,N,G,NI,NJ,NK,TMGRID,TNGRID,JCT,ICT 
          INTEGER :: MAXI,MAXJ,MINI,MINJ,ITEMP,JTEMP,IS,IE,JS,JE 
          REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: X,Y,Z
          CHARACTER(80) :: gridname
!
!         
!         Figure out min starting index and max ending indes
          ALLOCATE( GSI(NLGRIDS),GEI(NLGRIDS) )
          ALLOCATE( GSJ(NLGRIDS),GEJ(NLGRIDS) )
          DO G = 1,NLGRIDS
             MAXI = 0
             MAXJ = 0
             MINI = 99999
             MINJ = 99999
             DO N = 1,NLCASE
                IF (GNUMS(N) == G) THEN 
                   IF (IEA(N) > MAXI) THEN
                      MAXI = IEA(N) 
                   END IF
                   IF (JEA(N) > MAXJ) THEN
                      MAXJ = JEA(N)
                   END IF
                   IF (ISA(N) < MINI) THEN
                      MINI = ISA(N)
                   END IF
                   IF (JSA(N) < MINJ) THEN
                      MINJ = JSA(N)
                   END IF
                END IF 
             END DO 
             GSI(G) = MINI
             GEI(G) = MAXI
             GSJ(G) = MINJ
             GEJ(G) = MAXJ
          END DO 
!          WRITE(*,*) 'The max indicies',GSI(1),GEI(1),GSJ(1),GEJ(1) 
!          
!         GO through and real all the grid files already know how many
!         distinct grids we are using on this proc
!
          DO G = 1,NLGRIDS 
!            Allocate based on min/max indicies 
             ITEMP = GEI(G) - GSI(G) + 2          
             JTEMP = GEJ(G) - GSJ(G) + 2          
             ALLOCATE(GRIDIN(G)%X(ITEMP,JTEMP),GRIDIN(G)%Y(ITEMP,JTEMP))
             ALLOCATE(GRIDIN(G)%VOL(ITEMP,JTEMP)) 
             gridname = lgnames(G)
             TMGRID   = MGRID(G)  ! Number of grids in this file
             TNGRID   = NGRID(G)  ! Grid number we are interested in
             OPEN(FILE=gridname,UNIT=(RK+NPROCS+10),STATUS='OLD',&
                  FORM='UNFORMATTED') 
             READ((RK+NPROCS+10)) 
             IF (TNGRID .eq. 1) THEN
                READ((RK+NPROCS+10)) NI,NJ,NK
             ELSE
                READ((RK+NPROCS+10)) (NI,NJ,NK,N=1,TMGRID)  
             END IF    
!             WRITE(*,8) NI,NJ,NK
! 8           FORMAT('Reading Grid file, dimensioned ',3I7)
!            Store dimensions of array in module 
             NIA(G) = NI
             NJA(G) = NJ 
!            Now we know the dimesions of this particular grid (only
!            storing X,Y i.e. 2-D)
             ALLOCATE( X(NI,NJ,NK),Y(NI,NJ,NK),Z(NI,NJ,NK) )  
             DO N = 1,TMGRID-1
                READ((RK+NPROCS+10))             
             END DO             
             READ((RK+NPROCS+10)) (((X(I,J,K),I=1,NI),J=1,NJ),K=1,NK),& 
     &                            (((Y(I,J,K),I=1,NI),J=1,NJ),K=1,NK),& 
     &                            (((Z(I,J,K),I=1,NI),J=1,NJ),K=1,NK)
             CLOSE((RK+NPROCS+10))        
!
!            Only store the minimum number of grid pts GSI,GSJ contain
!            the starting I and J indicies for each grid
!             
!            We are assuming that we are running 2D so Z --> Y           
             JCT = 1
             DO J = GSJ(G),GEJ(G)
                ICT = 1
                DO I = GSI(G),GEI(G)
                   GRIDIN(G)%X(ICT,JCT) = X(I,J,2) 
                   GRIDIN(G)%Y(ICT,JCT) = Z(I,J,2) ! slight transformation
                   ICT = ICT + 1 
                END DO 
                JCT = JCT + 1
             END DO 
             DEALLOCATE( X,Y,Z )
          END DO ! end grid loop 
!              
          END SUBROUTINE READGRIDS 

!           
! 
!
      END MODULE GRID_MOD


