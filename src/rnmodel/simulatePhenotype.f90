subroutine SimulatePhenotype(verbose, nAnim, nComp, TBV, means, &
   vars, nLox, phen, locations, farmInd, pedigree)
! note that nfix and ncomp are different. nComp refers to number of components
! in simulation: for intercept and slope ncomp = 2. However, nfix depends on the
! type of analysis: if a single trait is desired nfix = 1, if random regression
! is to be conducted nfix = 2
  use constants, only: alloc1I, alloc1D, alloc2D, variances, &
   STDOUT, ZERO, STDERR, ONE
  use rng_module, only: gnormal
  use user_type
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nAnim, nComp, nLox
  real(KINDR), dimension(nAnim, nComp), intent(in) :: TBV
  type(variances), intent(in) :: vars
  real(KINDR), dimension(nComp), intent(in) :: means
  real(KINDR), dimension(nLox*nAnim), intent(out) :: phen
  real(KINDR), dimension(nanim,nlox), intent(out) :: locations
  integer, dimension(nanim*nlox), intent(out) :: farmInd
  integer, dimension(nanim, 3), intent(in) :: pedigree

  external :: covariate, defineFarms, allocateInd

  ! temporary variables
  integer :: i, nobs, iunfarm
  real(KINDR), dimension(1:ncomp,1:ncomp) :: temp2
  real(KINDR), dimension(1:ncomp) :: tempr ! to hold mean of zero
  real(KINDR), allocatable, dimension(:,:) :: E!, A, PE
  
  real(KINDR), dimension(:,:), allocatable, save :: farmBounds
  
  
  if (pedigree(1,1) == 0) then
   ! do nothing
  end if

  nobs = nAnim * nLox
  ! allocating temporary arrays
  !call alloc2D(X,nLox*nAnim,nFix, "locations", "simulatePhenotype")

  if (.not.allocated(farmBounds)) then
   call alloc2D(farmBounds, userParam%nFarm, 2, "farmBounds", "simulatePhenotype")
   call defineFarms(userParam%interval, userParam%nFarm, &
      userParam%farmRange, farmBounds)
   ! TODO: ask the user if they want this file
   open(newUnit = iunFarm, file = 'farms.txt')
   do i = 1, userParam%nfarm
     write(iunFarm, *) farmBounds(i, 1), farmBounds(i, 2)
   end do
   close(iunFarm)
  end if

  farmInd(1:nobs) = 0
  locations(1:nanim, 1:nlox) = ZERO

  ! allocating individuals, as in assingining each individual to a
  ! location and get "location" (challenge value) for that individual

  call allocateInd(nAnim, nlox, nobs, userParam%nFarm, userParam%allocation, &
   farmBounds, userParam%farmRange, farmInd, locations)

  ! covariance structure for E (covariances are zero)
  temp2(1:nComp, 1:nComp) = ZERO
  do i = 1, nComp
     temp2(i,i) = vars%E(i)
  end do
  tempr(1:nComp) = ZERO

  call alloc2D(E, nobs,nComp, "e", "simulatePhenotype")
  call gnormal(tempr, temp2, nComp, nobs, E)
  ! todo: implementation for PE if really it is required

  ! location is used instead of X because X is to be used in analysis
  ! and varies in dimension depending on the type of analysis
  call covariate(nComp, nObs, nAnim, nLox, TBV, E, phen, locations, means)
  if (verbose) write(STDOUT, *) " number of locations used:", nlox

end subroutine SimulatePhenotype

!!!! ============================================================ !!!!
subroutine allocateInd(nAnim, nlox, nobs, nfarm, allocation, farmBounds,&
     farmRange, farmInd, locations, pedigree, nm, male)
  use constants
  use rng_module
  implicit none
  integer, intent(in) :: nAnim, nLox, nobs, nFarm
  integer, intent(in) :: allocation
  real(KINDR), dimension(1:nfarm, 2), intent(in) :: farmBounds
  real(KINDR) :: farmRange
  integer, dimension(1:nobs), intent(out) :: farmInd
  real(KINDR), dimension(nAnim, nLox), intent(out) :: locations
  integer, dimension(1:nAnim, 3), optional :: pedigree
  integer, optional, intent(in) :: nm
  integer, optional, intent(in) :: male(:)

  integer, dimension(:), allocatable :: temp1, temp2
  integer :: spf! each farm has spf sires
  integer :: isire
  integer :: i, j, k
  real(KINDR) :: val
  real(KINDR), dimension(nlox) :: temp

  select case (allocation)
  case(1) ! random
     ! each individual is assigned to (only) one random farm
     do i = 1, nAnim
        call random_number(val)
        do while(val.eq.ONE)
           call random_number(val)
        end do
        j = int(val * nfarm) + 1
        k = nlox * (i - 1) + 1
        farmInd(k:(nlox*i)) = j
        call random_number(temp)
        ! shifting and scaling to match interval
        temp(1:nlox) = temp(1:nlox) * farmRange + farmBounds(j, 1)
        locations(i, 1:nlox) = temp(1:nlox)
        ! random_number is in [0,1] so no need to shift by 0 and divide by 1
     end do
  case(2) ! some sort of clustering:
     !one farm/sire & spf sires/farm & farms and sires random
     ! i.e., nf farms (as indices 1...nf) are randomised and each farm in this
     ! index array (temp2) will contain spf sires. The mating (pedigree) is 
     ! already decided.
     ! requires pedigree
     if (.not.present(pedigree)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " pedigree is required for allocation = 2"
        stop 2
     elseif (.not.present(nm)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " number of sires is required for allocation = 2"
        stop 2
     elseif (.not.present(male)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " male array is required for allocation = 2"
        stop 2
     elseif (nm < nfarm) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " number of farms is more than number of sires"
        stop 2
     end if
     spf = int(nm / nfarm)
     if (spf * nfarm .ne. nm) then
        write(STDERR,'(a)') "Error:"
        write(STDERR, *) "This is something that should not happen"
        write(STDERR, *) "The number of sires is not divisible to nfarms"
        stop 2
     end if

     call alloc1I(temp1, nm, "temp1", "allocateInd")
     call alloc1I(temp2, nm, "temp2", "allocateInd")
     do i = 1, spf
        j = (i-1)*nfarm + 1
        k = i * nfarm
        temp1(j:k) = (/(isire, isire = 1, nfarm)/)
     end do
     call choice(temp1, nm, nm, nm, temp2, nm) ! temp2 contains farms
     ! could do smarter way, but would require 2 more arrays
     !
     ! the method below requires the pedigree to be sorted, as this is the
     ! case for my simulations, I don't need to make this general.
     ! for the first individual, take the sire
     i = 1
     isire = pedigree(i,2)
     ! find its index in male array
     do j = 1, nm
        if (isire == male(j)) exit
     end do
     ! assign the farm to that individual
     k = temp2(j)
     farmind(((i-1)*nlox+1):(i*nlox)) = k
     call random_number(temp)
     locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
     do i = 2, nAnim
        !for all other individual, get the sire
        isire = pedigree(i,2)
        ! if its another offspring of the same sire, just assign the same
        ! again, assuming pedigree is sorted, this check is easy
        if (isire == pedigree(i-1, 2)) then
           k = temp2(j)
           farmind(((i-1)*nlox + 1):(i*nlox)) = k
           call random_number(temp)
           locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
           cycle ! and go to the next individual
        end if
        ! otherwise, find the index and do the rest
        do j = 1, nm
           if (isire == male(j)) exit
        end do
        k = temp2(j)
        farmind(((i-1)*nlox + 1):(i*nlox)) = k
        call random_number(temp)
        locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
     end do
     deallocate(temp1)
     deallocate(temp2)
  case(3:) ! later
  end select
  
end subroutine allocateInd

!!!! ============================================================ !!!!
subroutine defineFarms(interval, nfarm, diameter, farms)
  use constants
  use quickSort
  implicit none
  integer, intent(in) :: nfarm
  real(KINDR), dimension(2), intent(in) :: interval
  real(KINDR), intent(in) :: diameter
  real(KINDR), dimension(nfarm, 2), intent(out) :: farms

  real(KINDR), dimension(nfarm) :: centres
  integer, dimension(nfarm) :: ind
  
  call random_number(centres)
  call sortrx(nfarm, centres, ind)
  centres(1:nfarm) = centres(ind(1:nfarm))
  centres(1:nfarm) = (centres(1:nfarm) - centres(1)) / &
       (centres(nfarm) - centres(1)) * &
       (interval(2) - interval(1) - diameter) +&
       interval(1) + diameter / 2
  farms(1:nfarm, 1) = centres(1:nfarm) - diameter / 2
  farms(1:nfarm, 2) = centres(1:nfarm) + diameter / 2
  ! arithmetic floating point messes up with first and last boundary
  ! although experiment showed this is not required most of the times
  farms(1      , 1) = interval(1)
  farms(nfarm  , 2) = interval(2)
    
end subroutine defineFarms
