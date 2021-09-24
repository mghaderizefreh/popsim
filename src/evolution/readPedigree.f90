subroutine readPedigree(pedigreeFile, nanim, pedigree, nfounders)
use constants
implicit none
integer, intent(in) :: nanim
character(len = 100), intent(in) :: pedigreeFile
integer, dimension(nanim, 3), intent(inout) :: pedigree
integer, intent(out) :: nfounders

integer :: iunped
integer :: s, d, m, i
nfounders = 0
open(newUnit = iunped, file = trim(pedigreefile))
do i = 1, nanim
   read(iunped, *) s, d, m
   if ((m .lt. 0) .or. (d .lt. 0) .or. (s .le. 0)) then
      write(STDERR, '(a)') "Error:"
      write(STDERR, '(1x,a,i5)') "invalid input at line ", i
      write(STDERR, *) "Exiting"
      stop 1
   elseif ((m * d .eq. 0) .and. (m + d .ne. 0)) then
      write(STDERR, '(a)') "Error:"
      write(STDERR, '(1x,a,i5)') "only one parent is known for individual ", s
      write(STDERR, *) "not allowed at the moment. Exiting"
      stop 2
   elseif ((m + d .eq. 0)) then
      nfounders = nfounders + 1
   end if
   pedigree(i,1) = s
   pedigree(i,2) = d
   pedigree(i,3) = m
end do
end subroutine readPedigree
