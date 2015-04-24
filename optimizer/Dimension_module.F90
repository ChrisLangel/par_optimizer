      MODULE DIMEN_MOD
!
!
      IMPLICIT NONE
!
!      include 'mpif.h'
      INTEGER, SAVE ::RK,NPROCS,TOTALP,NUMCASES,IDEAL,NLCASE 
      INTEGER,ALLOCATABLE, DIMENSION(:,:), SAVE :: PMAP ! may want x
      INTEGER,ALLOCATABLE,DIMENSION(:),SAVE :: GRDPTS,ITA,JTA, & 
     &                                         ISA,JSA,IEA,JEA
      INTEGER,ALLOCATABLE,DIMENSION(:),SAVE ::GNUMS,GSI,GEI,GSJ, &
     &                                        GEJ,GIT,GJT
      CHARACTER(120),ALLOCATABLE,DIMENSION(:), SAVE :: cnames
      ! Watch the size here
      CHARACTER(80), DIMENSION(:), ALLOCATABLE,SAVE :: lcnames
!
      END MODULE DIMEN_MOD

