subroutine makeGenotype(nanim, nChr, nSNP, ident, genome, SNPlist, &
     chr_nlocibefore, genotypeFileBIN, genotypeFileTXT,&
     saveGenotypeBIN, saveGenotypeTXT)
  use constants
  use iso_fortran_env, only: int8
  implicit none

  integer, intent(in) :: nanim, nChr, nSNP !(nSNP per chromosome)
  integer, dimension(nanim), intent(in) :: ident
  type(chromosome), dimension(nChr), intent(in) :: genome
  integer, dimension(nChr,nSNP), intent(in) :: SNPlist
  integer, dimension(nchr), intent(in) :: chr_nlocibefore
  logical, intent(in) :: saveGenotypeBIN, saveGenotypeTXT
  integer(1), dimension(:), allocatable :: linear
  integer(8) :: linearSize, lsM1, counter, row, col
  character(len = 100), intent(in) :: genotypeFileTXT, genotypeFileBIN
  character(len = 100) :: formato

  integer :: totalSNP
  integer :: i, j, k
  integer :: a1, a2
  integer :: totLoci, ibit1, iblck1, iloci, id, iChr
  integer(kind=int8), dimension(:,:), allocatable :: genotype
  integer, dimension(:), allocatable, save :: pruningSNP
  integer, dimension(2,2) :: genocode

  interface
   subroutine cwritebin(a, nrow,ncol) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_long) :: nrow, ncol
    integer(c_signed_char) :: a(0:(nrow-1),0:(ncol-1))
   end subroutine cwritebin
  end interface
!  character(len = 100) :: formato

  genocode = reshape((/ 1, 2, 2, 3 /), (/ 2, 2 /))
  !genocode(1,1) = 1
  !genocode(1,2) = 2
  !genocode(2,1) = 2
  !genocode(2,2) = 3

  if (ident(1) > 0) then
  end if
  !--------------------------------------------
  ! reading genome to create the genotype
  !--------------------------------------------
  totLoci = 0
  do iChr = 1, nChr
     totLoci = totLoci + genome(iChr)%nloci
  end do
  totalSNP = nSNP * nChr
  if (.not.allocated(pruningSNP)) then
     call alloc1I(pruningSNP, totLoci, "pruningSNP", "GetGmatrix")
     pruningSNP(1:totLoci) = 0
     do iChr = 1, nchr
        pruningSNP(chr_nlociBefore(iChr) + SNPlist(iChr, :)) = 1
     end do
  end if
  allocate(genotype(nanim,totalSNP))
!  call alloc2I(genotype, nanim, totalSNP, "genotype", "GetGmatrix")

  do id = 1, nanim
     k = 0
     do iChr = 1, nChr

        i = 1
        iblck1 = 0 !block for loci i
        ibit1 = 0  ! bit for loci i
        do iloci = 1, genome(iChr)%nloci
           ibit1 = ibit1 - 1 ! position of the new loci i to be checked

           if (ibit1 < 0) then
              iblck1 = iblck1 + 1 ! a new block
              ibit1 = NBITS - 1 ! start from the 1st bit
           end if

           if (pruningSNP(chr_nlocibefore(iChr) + iloci) .eq. 0) cycle

           k = k + 1
           a1 = 1
           a2 = 1
           if (btest(genome(iChr)%genotypes(id, 1, iblck1), ibit1)) a1 = 2
           if (btest(genome(iChr)%genotypes(id, 2, iblck1), ibit1)) a2 = 2
           genotype(id, k) = genocode(a1, a2)
        end do
     end do
  end do
  !====================================
  if (saveGenotypeTXT) then
   open(1, file = trim(genotypeFileTXT))
   write(1, '(i0)') nanim
   write(1, '(i0)') totalSNP
   write(formato, '(a,i10,a)') "(i0,1x,",totalSNP-1,"(i1,1x),i1)"
   do i = 1, nanim
      write(1, formato) i, (genotype(i,j), j = 1, totalSNP)
   end do
   close(1)
  end if
  if (saveGenotypeBIN) then
   row = nanim
   col = totalSNP
   linearSize = row * col
   lsM1 = linearSize - 1
   allocate(linear(0:lsM1))
   do i = 0, (row - 1)
    do j = 0, (col - 1)
      counter = j
      counter = counter + i * col
      linear(counter) = genotype(i+1,j+1)
    end do
   end do
   call cwritebin(linear, row, col)
  end if
  
  deallocate(genotype)

end subroutine makeGenotype
