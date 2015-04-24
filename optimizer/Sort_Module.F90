! --------------------------------------------------------------------
! MODULE  Sorting:
!    This module can sort a set of numbers.  The method used is
! usually referred to as "selection" method.
! --------------------------------------------------------------------

MODULE  SORT_MOD
   USE DIMEN_MOD
   IMPLICIT  NONE
   PRIVATE   :: FindMinimum, Swap

CONTAINS

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)          ! assume the first is the min
      Location = Start             ! record its position
      DO i = Start+1, End          ! start with next elements
         IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
            Minimum  = x(i)        !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1             ! except for the last
         Location = FindMinimum(x, i, Size)  ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort


! --------------------------------------------------------------------
! SUBROUTINE  Sortdub():
!    This subroutine receives an array x() and sorts it into ascending
! order. Takes second array and swaps those elements in that order 
! --------------------------------------------------------------------

   SUBROUTINE  Sortdub(x, y, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(INOUT) :: x,y
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1             ! except for the last
         Location = FindMinimum(x, i, Size)  ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
         CALL  Swap(y(i), y(Location))  ! swap second array in same order 
      END DO
   END SUBROUTINE  Sortdub

! --------------------------------------------------------------------
! SUBROUTINE  Checkpop():
!   Subroutines tries to determine if it makes sense to yank a case 
!   off one processor and put on another
! --------------------------------------------------------------------

   SUBROUTINE Checkpop(IDIFF,TEMPA,PMAP,NCASES,I,J)  
      IMPLICIT  NONE
      INTEGER, INTENT(IN)  :: I,J
      INTEGER, DIMENSION(NPROCS), INTENT(INOUT) :: IDIFF,TEMPA,NCASES
      INTEGER, DIMENSION(NPROCS,NUMCASES), INTENT(INOUT) :: PMAP
      INTEGER, DIMENSION(NUMCASES) :: CHECK
      INTEGER :: N,P,M,PROC,LOC,BEST,IND,IDIS,NDIS,T1,T2,IND1
!
!
      BEST = 99999
      IF (IDIFF(I) < 0) THEN  ! see if we can pull one off
         DO P = 1,NPROCS      ! P is where we are placing the point
            IF (I .ne. P) THEN 
               CHECK = 0            ! zero out check array             
               DO N = 1,NCASES(TEMPA(I)) 
                  CHECK(N) = GRDPTS(PMAP(TEMPA(I),N)) + IDIFF(I) 
                  IF (ABS(CHECK(N)) + &
   &                 ABS(IDIFF(P) - GRDPTS(PMAP(TEMPA(I),N))) < BEST) THEN 
                     BEST = ABS(CHECK(N)) + & 
   &                 ABS(IDIFF(P) - GRDPTS(PMAP(TEMPA(I),N)))
                     IND  = N
                     PROC = P 
                  END IF
               END DO
            END IF  
         END DO 
      END IF 
!     See if placing this case on a different proc will balance the load better
      IDIS = ABS(IDIFF(I)) + ABS(IDIFF(J)) ! Initial # of unbalanced points
      NDIS = BEST 
      IF (NDIS <= IDIS) THEN 
         PMAP(TEMPA(PROC),NCASES(TEMPA(PROC))+1) = PMAP(TEMPA(I),IND) 
         DO N = IND,NUMCASES-1
            PMAP(TEMPA(I),IND) = PMAP(TEMPA(I),IND+1)
            IND = IND + 1              
         END DO 
         IDIFF(I)    = IDIFF(I) + GRDPTS(PMAP(TEMPA(I),IND))
         IDIFF(PROC) = IDIFF(PROC) - GRDPTS(PMAP(TEMPA(I),IND))
         NCASES(TEMPA(I)) = NCASES(TEMPA(I)) - 1
         NCASES(TEMPA(PROC)) = NCASES(TEMPA(PROC)) + 1
      END IF
   END SUBROUTINE Checkpop
      
     
!   Now look at swaps 
! --------------------------------------------------------------------
! SUBROUTINE  Checkswap():
!   Subroutines tries to determine if it makes sense to swap two cases
!   
! --------------------------------------------------------------------

   SUBROUTINE Checkswap(IDIFF,TEMPA,PMAP,NCASES,I,J)  
      IMPLICIT  NONE
      INTEGER, INTENT(IN)  :: I,J
      INTEGER, DIMENSION(NPROCS), INTENT(INOUT) :: IDIFF,TEMPA,NCASES
      INTEGER, DIMENSION(NPROCS,NUMCASES), INTENT(INOUT) :: PMAP
      INTEGER, DIMENSION(NUMCASES) :: CHECK
      INTEGER :: N,P,M,PROC,LOC,BEST,IND,IDIS,NDIS,T1,T2,IND1

      BEST = 99999
      DO P = 1,NPROCS
         IF (I .ne. P) THEN
            DO N = 1,NCASES(I)
                DO M = 1,NCASES(P)      
                   T1 = IDIFF(I) + GRDPTS(PMAP(I,N)) - &
    &                              GRDPTS(PMAP(P,M))              
                   T2 = IDIFF(P) - GRDPTS(PMAP(I,N)) + &
    &                              GRDPTS(PMAP(P,M))              
                   IF (ABS(T1) + ABS(T2) < BEST) THEN
                      BEST = ABS(T1) + ABS(T2) 
                      IND  = N
                      IND1 = M
                      PROC = P 
                   END IF  
               END DO 
            END DO
          END IF
      END DO 
      IDIS = ABS(IDIFF(I)) + ABS(IDIFF(PROC))  
!      WRITE(*,*) 'Hello',BEST,IDIS,I,PROC,GRDPTS(PMAP(I,IND)),&
!    &                                  GRDPTS(PMAP(PROC,IND1)) 
!      BEST = 9999
      IF (BEST < IDIS) THEN 
         T1               = PMAP(PROC,IND1)  
         PMAP(PROC,IND1)  = PMAP(I,IND)
         PMAP(I,IND)      = T1 
         T1 = IDIFF(I)    + GRDPTS(PMAP(I,IND)) - &
    &                       GRDPTS(PMAP(PROC,IND1))              
         T2 = IDIFF(PROC) - GRDPTS(PMAP(I,IND)) + &
    &                       GRDPTS(PMAP(PROC,IND1))              
          
         IDIFF(I)    = T1
         IDIFF(PROC) = T2
      END IF 
   END SUBROUTINE Checkswap 


END MODULE  SORT_MOD

