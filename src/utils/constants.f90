module constants
  implicit none

  ! Precision
  integer,        parameter :: KINDR = KIND(0d0)

  ! integer representation of alleles
  integer, parameter :: NBITS = 32

  !! derived types
  ! can make jagged arrays
  type JArr
     real(KINDR), dimension(:), allocatable :: array
  end type JArr

  type chromosome
     integer :: nloci, nblock
     real(KINDR):: chrL
     integer, dimension(:,:,:), pointer :: genotypes ! nanim x 2 x nblock
     ! if the loci are not equidist., then the position (no matter if Morgan or cM).
     real(KINDR), dimension(:), pointer :: positions ! nloci
  end type chromosome

  type QTL_Array
     integer :: nQTL ! number of QTL (on one chromosome)
     integer :: nComp ! number of traits being affected
     integer, dimension(:,:), allocatable :: indices ! nchar x nQTL
     real(KINDR), dimension(:,:,:), allocatable :: values ! nChr x nQTL x nComp
  end type QTL_Array

  type variances
   real(KINDR), allocatable, dimension(:) :: A ! genetic part
   real(KINDR), allocatable, dimension(:) :: E ! environemntal part
   real(KINDR), allocatable, dimension(:) :: PE ! permanent environment
   real(KINDR), allocatable, dimension(:,:) :: corr ! (genetic) correlation
   real(KINDR), allocatable, dimension(:,:) :: cov ! (genetic) covariance
  end type variances
  !   Handles
  integer, parameter :: STDIN  = 5
  integer, parameter :: STDOUT = 6
  integer, parameter :: STDERR = 0

  ! some other constants
  real(KINDR),    parameter :: ZERO   =  0_KINDR
  real(KINDR),    parameter :: ONE    =  1_KINDR
  real(KINDR),    parameter :: TWO    =  2_KINDR
  real(KINDR),    parameter :: FOURTH = .25_KINDR
  real(KINDR),    parameter :: HALF   =  .5_KINDR
  real(KINDR),    parameter :: THIRD  =  1._KINDR / 3._KINDR
  real(KINDR),    parameter :: PI     = 3.14159265358979323846_KINDR
  real(KINDR),    parameter :: SQRTPI = 1.7724538509055159_KINDR
  complex(KINDR), parameter :: I_C    = CMPLX(ZERO,    ONE,  KINDR)

  private :: alloc1D_L, alloc1D_I, alloc1D_D, alloc1D_J, alloc1D_Chr, &
    alloc1D_D_p
  private :: alloc2D_I, alloc2D_D, alloc2D_I8
  private :: alloc3D_I_p, alloc3D_D
  public :: alloc1D, alloc2D, alloc3D
  private :: HandleErr

  interface alloc1D  ! one dimenisonal arrays
    module procedure alloc1D_L ! logical
    module procedure alloc1D_I ! integer
    module procedure alloc1D_D ! real(double)
    module procedure alloc1D_D_p!real(double) pointer
    module procedure alloc1D_J ! jagged array
    module procedure alloc1D_Chr! chromosome
  end interface alloc1D

  interface alloc2D ! two dimensional (2D) array
    module procedure alloc2D_I ! integer
    module procedure alloc2D_I8! integer(kind8)
    module procedure alloc2D_D ! real(double)
  end interface alloc2D

  interface alloc3D ! 3D array or pointer
    module procedure alloc3D_I_p ! integer pointer
    module procedure alloc3D_D ! real(double)
  end interface alloc3D

contains
  subroutine alloc1D_L(array, len, name, from)
    implicit none
    logical, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D_L
  subroutine alloc1D_I(array, len, name, from)
    implicit none
    integer, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D_I

  subroutine alloc1D_D(array, len, name, from)
    implicit none
    real(KINDR), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D_D

  subroutine alloc1D_D_p(array, len, name, from)
    implicit none
    real(KINDR), dimension(:), pointer, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
  end subroutine alloc1D_D_p

  subroutine alloc1D_J(array, len, name, from)
    implicit none
    type(Jarr), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D_J

  subroutine alloc1D_Chr(array, len, name, from)
    implicit none
    type(chromosome), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D_Chr

  subroutine alloc2D_I(array, len1, len2, name, from)
    implicit none
    integer, dimension(:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc2D_I

  subroutine alloc2D_I8(array, len1, len2, name, from)
    use iso_fortran_env, only: int8
    implicit none
    integer(kind=int8), dimension(:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc2D_I8

  subroutine alloc2D_D(array, len1, len2, name, from)
    implicit none
    real(KINDR), dimension(:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2), stat = i)
       if (i .ne. 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc2D_D

  subroutine alloc3D_I_p(array, len1, len2, len3, name, from)
    implicit none
    integer, dimension(:,:,:), pointer, intent(inout) :: array
    integer, intent(in) :: len1, len2, len3
    character(len=*), intent(in) :: name, from
    integer :: i
    allocate(array(len1, len2, len3), stat = i)
    if (i .ne. 0) call HandleErr(i, name, from)
  end subroutine alloc3D_I_p
  subroutine alloc3D_D(array, len1, len2, len3, name, from)
    implicit none
    real(KINDR), dimension(:,:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2, len3
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2, len3), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc3D_D
  
  subroutine HandleErr(i, name, from)
    implicit none
    character(len=*), intent(in) :: name, from
    integer, intent(in) :: i
    if (i == -1) then
       write(STDERR, '(A)') "Error"
       write(STDERR,*) "array", trim(name), " in ", "is already allocated"
       stop 2
       if (i == 1) then
          write(STDERR,'(A)') "Error in system routine:"
          write(STDERR,*) "cannot allocate ", trim(name), " in ", trim(from)
          stop 2
       elseif (i == 2) then
          write(STDERR, '(A)') "Error (invalid data):"
          write(STDERR,*) "cannot allocate", trim(name), " in ", trim(from)
          stop 2
       elseif (i == 3) then
          write(STDERR, '(A)') "Error (invalid data and sys routine):"
          write(STDERR,*) "cannot allocate", trim(name), " in ", trim(from)
          stop 2
       else
          write(STDERR, '(A)') "Error:"
          write(STDERR, *) "cannot allocate", trim(name), " in ", trim(from),&
               "for unknown reason"
          stop 2
       end if
    end if
  end subroutine HandleErr
     
  
end module constants

