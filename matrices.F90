!----------------------------------------------------------------------
!
!   module name        - matrices
!
!----------------------------------------------------------------------
!
!   latest revision    - May 2
!
!   purpose            - store element matrices for the first layer
!                        elements
!
!----------------------------------------------------------------------
!
module matrices
!
  use parametersDPG
  implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!
! max # of elements in the first layer
  integer, parameter :: MAXNRFL=4
!
! # of stored elements in the first layer
  integer :: NRFL
!!$OMP THREADPRIVATE (NRFL)
!
! order of elements
  integer, parameter :: MYP=6
!
! matrix dimensions
  integer, parameter :: MYE = 3*MYP*(MYP+1)**2*2
  integer, parameter :: MYQ = MYP**3*6
!
! stiffness matrices to store
  VTYPE, dimension(MYE,MYE,MAXNRFL) :: ZFL_EE
!!$OMP THREADPRIVATE (ZFL_EE)
  VTYPE, dimension(MYE,MYQ,MAXNRFL) :: ZFL_EQ
!!$OMP THREADPRIVATE (ZFL_EQ)
  VTYPE, dimension(MYQ,MYQ,MAXNRFL) :: ZFL_QQ
!!$OMP THREADPRIVATE (ZFL_QQ)
  VTYPE, dimension(MYE,MAXNRFL) :: ZLOADFL_E
  VTYPE, dimension(MYQ,MAXNRFL) :: ZLOADFL_Q
!
! xy coordinates of the first vertex node
  real*8, dimension(2,MAXNRFL) :: XYVERT
!!$OMP THREADPRIVATE (XYVERT)
!
end module matrices



