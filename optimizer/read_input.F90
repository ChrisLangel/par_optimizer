SUBROUTINE READ_INPUT

USE input
USE OPT_MOD
USE PARAM_MOD

IMPLICIT NONE
!
TYPE NODE
   TYPE(NODE), POINTER :: next
   REAL(KIND=8)        :: VAL
END TYPE NODE 

TYPE(NODE), POINTER :: PARAMS,FIRST 

LOGICAL :: eof,eol
CHARACTER(LEN=16) :: w
DOUBLE PRECISION :: T1,T2
INTEGER :: i,spaces,nump

call input_options(echo_lines=.false.,error_flag=1)

eof = .false.
open(file='param_inp',unit=2)
do while (eof .eqv. .false.) ! read until end of file 
call read_line(eof,2)
! Read first word to determine what parameter we are setting
! will automatically skip commented lines 
call readu(w) 
if (w .ne. " ") then 
   select case(w)  ! figure out what parameter is what  
   case("OPTYPE")
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w)   ! do loop allows varrying spacings (e.g. P=1, P =1, P = 1) 
       end do   
       otype = w
   case("STEPFRAC")
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) STEPFRAC ! convert to double precision 
   case("MAXITS") 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) MAXITS 
   case("TOL") 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) TOL
   case("DIFFSTEP") 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) DIFFSTEP
   case("INTNCN") 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) INTNCN 
   case("PTNCON") 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       read (w,*) PTNCON
   case("PMVEC")  ! a bit tricky to make it accept any length vector 
       call readu(w)     
       do while (w .eq. " ") 
       call readu(w) 
       end do
       allocate( FIRST ) 
       read (w,*) T1
       FIRST%VAL =  T1
       PARAMS    => FIRST 
       nump   = 1
       spaces = 0 
       eol    = .true. 
       do while (eol) 
          call readu(w) 
          if (w .ne. " ") then
            read (w,*) T1
            allocate( PARAMS%next )
            PARAMS%next%VAL =  T1
            PARAMS          => PARAMS%next 
            nump = nump + 1      
            spaces = 0  
          else
            spaces = spaces + 1
          end if
          if (spaces > 3) then ! say if there is more than 3 spaces 
            eol = .false.      ! we have hit the end of the line
          end if
       end do 
       nullify( PARAMS%next ) 
   end select 
end if 
end do 


NCON = nump 
allocate( PMVEC(NCON) ) 
PARAMS => FIRST
do i = 1,NCON
   PMVEC(i)  =  PARAMS%VAL
   PARAMS    => PARAMS%next 
end do 

! allocate gradient vectors 
allocate( GRADVEC(NCON) )
allocate( STEPDIR(NCON),PMVECK(NCON),GRADVECK(NCON) )
allocate( Hk(NCON,NCON) )


END SUBROUTINE READ_INPUT 
