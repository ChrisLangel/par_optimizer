      SUBROUTINE SORTCASES
!
      USE DIMEN_MOD
      USE SORT_MOD 
!
      IMPLICIT NONE
!     
      TYPE NODE
         TYPE(NODE), POINTER :: next
         INTEGER             :: IND
         INTEGER             :: PNUM  
      END TYPE NODE 
!      
!    
      TYPE(NODE), POINTER   :: PROCMAP,FIRST
! 
      INTEGER :: I,J,N,P,ICT,PCT,IND,TEMP,HALF,NUMLEFT,DUBS
      INTEGER :: TPROC,SCORE,TSCORE,FSCORE 
      INTEGER,DIMENSION(NPROCS) :: NPTS,TEMPA,IDIFF,NCASES,TIDIFF,FIDIFF
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: TPMAP,FPMAP
      LOGICAL :: SORTED
      CHARACTER(120) :: tempstring,filename 
!
! 
      ALLOCATE ( PMAP(NPROCS,NUMCASES) )
      ALLOCATE ( TPMAP(NPROCS,NUMCASES))
      ALLOCATE ( FPMAP(NPROCS,NUMCASES))
! 
!     Sort the cases by how many points per case (bubble sort)
      SORTED = .FALSE. 
      DO WHILE (SORTED .eqv. .FALSE.)
      SORTED = .TRUE. 
         DO I = 1,NUMCASES-1
            IF ( GRDPTS(I+1) < GRDPTS(I) ) THEN 
               TEMP        = GRDPTS(I+1)
               tempstring  = cnames(I+1)   
               GRDPTS(I+1) = GRDPTS(I)
               GRDPTS(I)   = TEMP
               cnames(I+1) = cnames(I) 
               cnames(I)   = tempstring
               SORTED      = .FALSE. 
            END IF 
         END DO 
      END DO
!
      DO I = 1,NPROCS
         TEMPA(I) = I
      ENDDO 
!
      IDEAL = TOTALP/NPROCS
      I       = 1
      J       = NUMCASES
      NUMLEFT = NUMCASES
      ICT     = 1
      PCT     = 1 
      NPTS    = 0 
      PMAP    = 0  
      ! We should only do this loop if there is more than one proc
      IF (NPROCS > 1) THEN
         DO WHILE (J > I .and. NUMLEFT/2 >= (NPROCS-PCT+1)) ! make sure we can load two on each remaining proc
            PMAP(PCT,ICT)   = I
            PMAP(PCT,ICT+1) = J
            NPTS(PCT) = NPTS(PCT) + GRDPTS(I) + GRDPTS(J) 
            PCT = PCT + 1 
            IF (PCT > NPROCS) THEN 
               PCT = 1  
               ICT = ICT + 2
            END IF   
            J = J - 1 
            I = I + 1
            NUMLEFT = NUMLEFT - 2
         END DO
!        
!        Sort by number of points on each proc        
         CALL Sortdub(NPTS,TEMPA,NPROCS)
!        For remaining cases distribute based on number current number of points on proc
!        TEMPA is an aray of proc #'s ordered by how many points on that proc
         IND   = 1
         TPROC = NPROCS
         DUBS  = NUMLEFT - NPROCS ! # of procs that will have two cases loaded
         DO WHILE (NUMLEFT > 0) 
            !WRITE(*,*) J
            IF (NUMLEFT > TPROC) THEN ! load two on proc with fewest
               PMAP(TEMPA(IND),ICT)   = I
               PMAP(TEMPA(IND),ICT+1) = J-NPROCS+DUBS 
               I       = I + 1
               DUBS    = DUBS  - 1
               TPROC   = TPROC - 1
               IND     = IND + 1
               NUMLEFT = NUMLEFT - 2
            ELSE
               PMAP(TEMPA(IND),ICT)   = J 
               J       = J - 1
               IND     = IND + 1
               NUMLEFT = NUMLEFT - 1 
            END IF
         END DO       
!           
!        Check how many points and cases are on each proc
!        Also check how close to ideal we have gotten
!          
         CALL GETSTATS(NPTS,NCASES,IDIFF)

         SCORE  = 0
         DO I = 1,NPROCS
            SCORE = SCORE + ABS(IDIFF(I))  
!            WRITE(*,*) IDIFF(I),NPTS(I),PMAP(I,1),PMAP(I,2),PMAP(I,3)
         END DO
         FSCORE = SCORE 
!         WRITE(*,*) FSCORE   
         FPMAP  = PMAP 
         TPMAP  = PMAP
         TIDIFF = IDIFF
!
!        This bit will do a second pass and try and even out the load as best as possible
         DO P = 1,1 
            DO I = 1,NPROCS
               TEMPA(I) = I
            ENDDO 
            CALL Sortdub(IDIFF,TEMPA,NPROCS)
            DO N = 1,NPROCS
               CALL Checkpop(IDIFF,TEMPA,PMAP,NCASES,N,J)
               TSCORE = 0
               DO I = 1,NPROCS
                  TSCORE = TSCORE + ABS(IDIFF(I))
               END DO 
               IF (TSCORE <= SCORE) THEN ! Make sure it is actually better
                  SCORE  = TSCORE
                  TIDIFF = IDIFF
                  TPMAP  = PMAP
               ELSE 
                  PMAP   = TPMAP
                  IDIFF  = TIDIFF
                  CALL GETSTATS(NPTS,NCASES,IDIFF)  
                  DO I = 1,NPROCS
                     TEMPA(I) = I
                  ENDDO 
                  CALL Sortdub(IDIFF,TEMPA,NPROCS)
               END IF
            END DO  
         END DO         
         CALL GETSTATS(NPTS,NCASES,IDIFF)

         DO I = 1,NPROCS
            TEMPA(I) = I
         ENDDO
         N = 1 
         CALL Sortdub(IDIFF,TEMPA,NPROCS)
         CALL Checkpop(IDIFF,TEMPA,PMAP,NCASES,N,J) 
         TSCORE = 0
         DO I = 1,NPROCS
             TSCORE = TSCORE + ABS(IDIFF(I))
         END DO
         IF (TSCORE <= SCORE) THEN ! Make sure it is actually better
            SCORE  = TSCORE
            TIDIFF = IDIFF
            TPMAP  = PMAP
         ELSE 
            PMAP   = TPMAP
            IDIFF  = TIDIFF
         END IF
         CALL GETSTATS(NPTS,NCASES,IDIFF)
!        ELSE !We are on single proc
         CALL Checkswap(IDIFF,TEMPA,PMAP,NCASES,N,J) 
         CALL GETSTATS(NPTS,NCASES,IDIFF)
         TSCORE = 0 
         DO I = 1,NPROCS
            TSCORE = TSCORE + ABS(IDIFF(I))
         END DO
!         WRITE(*,*) 'This is the temp score', TSCORE 
         IF (TSCORE > FSCORE) THEN
            PMAP = FPMAP    
         END IF
         CALL GETSTATS(NPTS,NCASES,IDIFF)
      ELSE ! we are just on one proc 
         NCASES(1) = NUMCASES
         DO N = 1,NUMCASES
            PMAP(1,N) = N 
         END DO         
      END IF ! 
!
    
      DO I = 1,NPROCS
         WRITE(*,*) IDIFF(I),NPTS(I),PMAP(I,1),PMAP(I,2),PMAP(I,3),&
     &         GRDPTS(PMAP(I,1)),GRDPTS(PMAP(I,2)),GRDPTS(PMAP(I,3)) 
      END DO    

!     Write out files with info about cases on blocks
      DO I = 1,NPROCS
         WRITE (filename, '("Pcases", I2.2, ".txt")' )  I
         OPEN (FILE=filename,UNIT=I+10,STATUS='REPLACE',FORM='FORMATTED')
         WRITE((I+10),10) NCASES(I) 
         DO J = 1,NCASES(I)
            WRITE((I+10),11) cnames(PMAP(I,J)) 
         END DO 
         CLOSE((I+10))
      END DO 
 10   FORMAT(I5)
 11   FORMAT(A)
      OPEN(999,FILE='donesplit',STATUS='REPLACE')


      END SUBROUTINE
!
!     Subroutine that returns the number of points/number of
!     cases on each proc   
! 
      SUBROUTINE GETSTATS(NPTS,NCASES,IDIFF)
      USE DIMEN_MOD
      IMPLICIT NONE
      INTEGER, DIMENSION(NPROCS),INTENT(INOUT) :: NPTS,NCASES,IDIFF
      INTEGER                             :: ICT,IND,I            

         NPTS   = 0      ! Zero out counters
         NCASES = 0
         DO I = 1,NPROCS
            ICT = 1
            IND = 1
            DO WHILE (IND > 0)
               NPTS(I)   = NPTS(I) + GRDPTS(PMAP(I,ICT))
               NCASES(I) = NCASES(I) + 1 
               ICT       = ICT + 1
               IND       = PMAP(I,ICT) ! see if next entry is zero
            END DO
            IDIFF(I) = IDEAL - NPTS(I)     
         END DO  


       END SUBROUTINE GETSTATS
