      MODULE PARAM_MOD
!
!
      IMPLICIT NONE
!
!      include 'mpif.h'
      INTEGER, SAVE :: NCON,OPTITER 
      REAL(KIND=8)  :: LFVAL,TFVAL 
!     The vector of coefficients we are optimizing    
      REAL(KIND=8),DIMENSION(:),ALLOCATABLE, SAVE :: PMVEC
!
      END MODULE PARAM_MOD

