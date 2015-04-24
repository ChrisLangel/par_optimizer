!
!
      SUBROUTINE GETGRIDCORDS
!
      USE GRID_MOD
      USE DIMEN_MOD
!
      IMPLICIT NONE
!
      INTEGER             :: IO,TEMP,N,I,J,MG,NG
      INTEGER             :: IT,JT,LT,IS,JS,IE,JE
      INTEGER, DIMENSION(:), ALLOCATABLE :: TMGRID,TNGRID
      CHARACTER(80)       :: filename,tempname
      CHARACTER(80)       :: gridname ! may need to change size of this
      CHARACTER(80), DIMENSION(:), ALLOCATABLE :: tgnames,tgnamesU
      LOGICAL             :: EXISTS 
!
      IO = 1
!
!     Make sure the root processor has fininshed writing the connectivity
!     files.  If not just loop until it has.
!  
 100  CONTINUE     
      INQUIRE(FILE='donesplit',EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         GOTO 100
      ENDIF
!
!     Open case file that is associated with this processor  
      WRITE (filename, '("Pcases", I2.2, ".txt")' )  RK+1 
      OPEN (FILE=filename,UNIT=(RK+100),STATUS='OLD', &
     &       FORM='FORMATTED')
      REWIND(RK+100)  
!
!     First line tells us how many cases 
      READ((RK+100),*) TEMP
      NLCASE = TEMP
!     
!     Allocate array with all the case names and local grid array
      ALLOCATE( lcnames(TEMP) ) 
      ALLOCATE( tgnames(TEMP),tgnamesU(TEMP) ) 
      ALLOCATE( TMGRID(TEMP),TNGRID(TEMP) ) 
!     Arrays that store the indicies that are used for each case
      ALLOCATE( GIT(TEMP),GJT(TEMP),ISA(TEMP) )
      ALLOCATE( JSA(TEMP),IEA(TEMP),JEA(TEMP) ) 
!     Array that stores grid number of each case
      ALLOCATE( GNUMS(TEMP) ) 
! 
!     Go through cases on this proc and determine what grids will be
!     used for each case       
      DO N = 1, NLCASE
         READ((RK+100),10,IOSTAT=IO) tempname   
         OPEN(FILE=tempname,UNIT=(RK+NPROCS+10),STATUS='OLD',&
     &        ACTION='READ',FORM='UNFORMATTED') 
         REWIND((RK+NPROCS+10)) 
         READ((RK+NPROCS+10)) gridname
         lcnames(N) = tempname
         tgnames(N) = gridname  
         READ((RK+NPROCS+10)) MG,NG
         READ((RK+NPROCS+10)) IT,JT,LT,IS,IE,JS,JE
         TMGRID(N) = MG
         TNGRID(N) = NG
         GIT(N)    = IT
         GJT(N)    = JT
         ISA(N)    = IS
         IEA(N)    = IE
         JSA(N)    = JS
         JEA(N)    = JE
!         WRITE(*,10) tempname,gridname
         CLOSE((RK+NPROCS+10))
      END DO 
 10   FORMAT(A) 
 11   FORMAT(I5) 
      ! Ditch the actual local case file 
!      CLOSE((RK+100),STATUS="DELETE") 
!     Create a temporary array that is not changed in next loop       
      tgnamesU = tgnames
!
!     Look for repeat grids 
      DO I = 1,NLCASE
         DO J = 1,NLCASE
            IF (I .ne. J) THEN
               IF (tgnames(I) .eq. tgnames(J) ) THEN
                  tgnames(I) = 'same' 
               END IF 
            END IF
         END DO
      END DO 
!
!     Create an array of unique grids 
      N = 0
      DO I =1,NLCASE ! count number of unique grids 
         IF (tgnames(I) .ne. 'same') THEN 
            N = N + 1       
         END IF
      END DO
      NLGRIDS = N    ! store in module and allocate arrays
      ALLOCATE( lgnames(NLGRIDS) )   
      ALLOCATE( MGRID(NLGRIDS),NGRID(NLGRIDS) )
      N = 1 
      DO I =1,NLCASE ! map to stored array grid name and multigrid info
         IF (tgnames(I) .ne. 'same') THEN 
            lgnames(N) = tgnames(I)
            NGRID(N)   = TNGRID(I)
            MGRID(N)   = TMGRID(I)    
            N = N + 1  
         END IF
      END DO
!     Get rid of temporary arrays
      DEALLOCATE( TMGRID,TNGRID )  
      DEALLOCATE( tgnames )
!
!     Go through and figure out which number grid each case is using
      DO I = 1,NLCASE
         DO J = 1,NLGRIDS 
            IF (tgnamesU(I) .eq. lgnames(J)) THEN
               GNUMS(I) = J
            END IF
         END DO
      END DO 
!   
!     Allocate the main grid array for this proc
      IF (.NOT. ALLOCATED(GRIDIN) ) THEN 
         ALLOCATE( GRIDIN(NLGRIDS) )
      ENDIF
      ALLOCATE( NIA(NLGRIDS),NJA(NLGRIDS) ) 
!
!     Call function that reads the grids and store coordinates in GRIDIN      
      CALL READGRIDS     
!
      END SUBROUTINE

