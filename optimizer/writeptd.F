C
C     Program that writes data in P3d format 
C     HEADER: Jtot,Ktot,Ltot,NumQvar 
C      
      PROGRAM WRITEQ
C
C
      INTEGER :: NJ,NK,NL,NQ,J,K,L,N
      INTEGER :: IT,JT,IS,JS,IE,JE
      REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Q
      CHARACTER(80) :: OFILE,IFILE,gname,testname
C       
      NJ = 3
      NK = 4
      NL = 5  
      NQ = 4 
      IT = 412
      JT = 160
      IS = 220
      IE = 250
      JS = 20
!      JE = 40


      WRITE(*,20,ADVANCE='NO') 'Enter name of grid file: '
      READ(*,20,END=100) gname
 20   FORMAT(A)

      WRITE(*,20,ADVANCE='NO') 'Enter name of test: '
      READ(*,20,END=100) testname


      WRITE(*,20,ADVANCE='NO') 'Enter end of J (> 21): '
      READ(*,*,END=100) JE

      OPEN(3,FILE=testname,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(3) gname
      WRITE(3) 1,1  
      WRITE(3) IT,JT,IS,IE,JS,JE
C
      ALLOCATE(Q(IT,JT,NQ))
      DO N=1,NQ
      DO K=1,JT
      DO J=1,IT 
         Q(J,K,N) = REAL(N)*REAL(K)
      ENDDO
      ENDDO
      ENDDO
      WRITE(3) Q
      CLOSE(3)      
100   CONTINUE       
      END      
