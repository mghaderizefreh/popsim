subroutine gnormal(mean, cov, dim, N, output, seed)
  ! this subroutine produces normal distribution using boxMuller transformation


  ! Inputs:
  !     `dim`    : dimension of the distribution (integer)
  !     `mean`   : mean of the distribution (1D real(KINDR) array of shape `dim`)
  !     `cov`    : covariance matrix (2D real(KINDR) array of shape `dim` x `dim`)
  !     `N`   : number of samples (integer)
  !     `output` : (2D real(KINDR) array of shape `dim` x `N`) output
  !     `nseed`  : seed for random stream (integer)
  ! matrix `cov`. `N` is the length of the array and `output` is a real(KINDR)
  ! array of shape `dim` x `output`.
  !
  ! written by Masoud Ghaderi Zefreh
  ! first revision : 1 December 2020
  
  use constants

  implicit none
  integer, intent(in) :: dim, N
  real(KINDR), dimension(1:dim), intent(in) :: mean
  real(KINDR), dimension(1:dim, 1:dim), intent(in) :: cov
  real(KINDR), dimension(1:N, 1:dim), intent(inout) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer :: dim2, i, j
  real(KINDR) :: sd, val
  real(KINDR), dimension(dim, dim) :: covcopy
  real(KINDR), dimension(:,:), allocatable :: uniform

  external :: dpotrf, dgemm, normdev

  covcopy = cov

  if (present(seed)) then
     call random_seed(put = seed)
  end if

  if (dim == 1) then
     sd = sqrt(cov(1,1))
     do i = 1, N
        call normdev(val)
        output(i,1) = mean(1) + sd * val
     end do
     return
  else

     dim2 = dim
     if (mod(dim, 2) .eq. 1) dim2 = dim2 + 1

     ! making uniform distribution
     call alloc2D(uniform, N, dim2, "uniform", "gnormal")
     call random_number(uniform)

     ! making independent pairs of normally distributed data
     do i = 1, dim, 2
        j = i + 1
        output(:,i) = sqrt(-2 * log(uniform(:,i))) * cos(2 * PI * uniform(:,j))
     end do
     do j = 2, dim, 2
        i = j - 1
        output(:,j) = sqrt(-2 * log(uniform(:,i))) * sin(2 * PI * uniform(:,j))
     end do

     if (dim .ne. dim2)then
        deallocate(uniform)
        call alloc2D(uniform, N,dim, "uniform", "gnormal")
     end if
     uniform = output

     ! cholesky factorisation
     call dpotrf('L',dim,covcopy,dim,i)
     if (i.gt.0) then
        write(STDERR, *) ' Cholesky factorisation failed!'
        write(STDERR, *) ' Covariance matrix is not positive definite'
        write(STDERR, *) ' Exiting...'
        stop 2
     elseif (i.lt.0) then
        write(STDERR, *) "some other shit happened"
        stop 2
     end if

     ! the cholesky factorisation routine ?potrf(with 'L') returns the lower
     ! part of the cholesky. The upper part must be zero.
     do i = 1, dim - 1
        do j = i + 1, dim
           covcopy(i,j) = ZERO
        end do
     end do

     ! the output is normal(nxd) = unifrm(nxd) * chol.Transpose
     call dgemm('N','T', N, dim, dim, ONE, uniform, N, covcopy, dim, &
          ZERO, output, N)

     forall( i = 1:N)
        output(i,1:DIM) = mean(1:DIM) + output(i,1:DIM)
     end forall
  end if
  deallocate(uniform)
end subroutine gnormal

subroutine normdev(val)
! making a random varaible form standard normal
! copy of Andrea's code (excluding sigma)
  use constants
  implicit none
  real(KINDR), intent(inout) :: val
  real(KINDR) :: r, b, x, y

  call random_number(x)
  call random_number(y)
  
  r = sqrt(-TWO * Log(x))
  b = y * 6.2831853_KINDR
  
  val = r * cos(b)
end subroutine normdev
