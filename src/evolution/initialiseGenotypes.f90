!==============================================================================================
! this subroutine initialise the genotype for the first generation of the historical generations
!
! there are tree options for setting genotypes (genstart parameter)
! genstart=1  =  all loci set to be fixed to 
! genstart=2  =  all loci are with frequ 0.5, in HWE and in linkage equilibrium (i.e. genotypes are sampled independently to be in HWE )
! genstart=3  =  genotypes are read from files
!
!if option 1 and 2 are used, all chromosomes have similar 
!

subroutine InitialiseGenotypes(verbose, nchr, nanim, genstart, nloci, nblock,&
     istore, nfounders, genome, maxloci, maxblock, ifail, chrL,&
     prefixfilename1, prefixfilename2)
  use constants
  use rng_module
  use quickSort
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: genstart, nanim, nchr, istore, nfounders
  character(len=*), optional, intent(in) :: prefixfilename1, prefixfilename2
  integer, intent(inout) :: nloci
  type(chromosome), dimension(nChr), intent(inout) :: genome
  real(KINDR), intent(in) :: chrL
  integer, intent(out) :: ifail, nblock
  integer :: iun, i, k, id, j, ichr, igam, a1, iblock, nlines
  real(KINDR) :: rand
  integer :: maxblock, maxloci
  character(len=256) :: filename
  integer, dimension(:), allocatable :: temp, animIndex, gamindex, binindex

  !nloci is input if genstart is 1 or 2
  if (verbose) write(STDOUT, '(a35, i1,a1)') " Initialising genotypes&
       &(genstart: ",genstart, ")"
  if (genstart .eq. 1 .or. genstart .eq. 2) then
     write(STDERR, '(a)') "ERROR:"
     write(STDERR, '(a)') "implementation needed for positions"
     stop 2
  end if
  !=====================================================
  ! IF GENSTART=1 or 2
  ! all chrmosmoe have the same number of SNP
  !
  ! the number of loci is passed as input in dummy variable <nloci>
  !=====================================================
  !
  ! This block of code initialise the array of genome TYPE so that
  ! each element has its 'genotypes' pointer array property
  ! allocated to nanim x 2 x nblock
  if (genstart .eq. 1 .or. genstart .eq. 2) then
     if(istore .eq. 1) then
        i = 32
        nblock = nloci/i
        if(mod(nloci,i) > 0) nblock = nblock + 1
     else
        nblock = nloci
     endif
     maxblock = nblock
     maxloci = nloci
     do ichr = 1, nchr
        genome(ichr)%nloci  = nloci
        genome(ichr)%nblock = nblock
        genome(ichr)%chrL = chrL
        call alloc3D(genome(ichr)%genotypes, nanim, 2, nblock, &
             "genome(ichr)%genotypes", "initialiseGenotypes")
     end do
  end if
  ifail = 0

  !=====================================================
  ! GENSTART=1
  ! ALL SNP in ALL CHR SET TO BE  FIXED 
  !
  !=====================================================
  if (genstart .eq. 1) then
     k = 0
     if (istore .eq. 1) then
        do j = 31, 0, - 1
           k = ibclr(k, j ) !ibclr(k,j) returns the value of k  with the bit at position j set to zero.
        end do
     else
        k = 1  !genotype are 1
     end if
     do ichr = 1, nchr
        genome(ichr)%genotypes(:,:,:) = k !all animalss are homozygous for all genes
     end do

     !=====================================================
     ! GENSTART=2
     !  SNP ARE SEGREGATING WITH FREQ =0.5 IN HWE AND IN LINKAGE EQUILIBRIUM
     !
     !=====================================================
  elseif (genstart .eq. 2) then
     ! if the initial genotype are segregating at freq 0.5 and in LE across SNPs
     ! when genotypes are store in bits
     if (istore .eq. 1) then
        do ichr = 1, nchr
           do iblock = 1, nblock
              do igam = 1, 2
                 do id = 1, nanim
                    a1 = k
                    do j = 31, 0, - 1
                       call random_number(rand)
                       !ibset(i,j) returns the value of i with 
                       !the bit at position j set to one.
                       if ( rand > HALF ) a1 = ibset( a1, j ) 
                    end do
                    genome(ichr)%genotypes(id, igam, iblock ) = a1
                 end do
              end do
           end do
        end do
     else
        !  when genotypes are store as integer
        a1 = 1
        do ichr = 1, nchr
           do iblock = 1, nblock
              do igam = 1, 2
                 do id = 1, nanim
                    call random_number(rand)
                    !change the allele (it is ok as frequency
                    ! is the same of both alleles)
                    if ( rand > 0.5 ) a1 = 3 - a1  
                    genome(ichr)%genotypes(id, igam, iblock ) = a1
                 end do
              end do
           end do
        end do
     end if

     !=====================================================
     ! GENSTART=3
     ! SNP GENOTYPE ARE READ FROM FILES
     ! CHROMOSOMES CAN HAVE DIFFERENT NUMBER OF SNP/LOCI
     !
     !=====================================================
  elseif (genstart .eq. 3) then

     maxblock = 0
     maxloci = 0
     chr: do ichr = 1, nChr
        genome(ichr)%chrL = chrL
        write(filename,'(a,i3.3)') trim(prefixfilename1), ichr
        if (verbose) write(STDOUT,'(a, a)') "  reading ", trim(filename)
        open(newunit = iun, file=filename, status = 'old')
        read(iun,*) i, nloci, nblock, k

        if (i .lt. nfounders) then
           ifail = 1
           write(STDERR,'(a, a)') &
                "Error. Genotype file does not have enough individuals: ",&
                trim(filename)
           write(STDERR,'(a, 2i8)') "  i nanim", i, nanim
           stop 2
        end if
        if (nloci > maxloci) maxloci = nloci
        if (nblock > maxblock) maxblock = nblock

        genome(ichr)%nloci  = nloci
        genome(ichr)%nblock = nblock
        call alloc3D(genome(ichr)%genotypes, nanim, 2, nblock, &
             "genome(ichr)%genotypes", "initialiseGenotypes")

        nLines = i * 2 ! number of lines
        i = 2 * nFounders
        if (.not.allocated(temp)) then
           call alloc1D(temp, nLines, "temp", "initialiseGenotypes")
           call alloc1D(gamIndex, i, "index", "initialiseGenotypes")
           call alloc1D(binIndex, i, "gamIndex","initialiseGenotypes")
           call alloc1D(animIndex, i, "gamIndex","initialiseGenotypes")
        end if

        temp = (/( id, id = 1, nLines )/)
        !gamindex = ind of lines to read
        call choice(temp, nLines, nLines, i, gamIndex, i)
        ! index of animals (1..(2xnFounders)) to be genotypes in random order
        call choice(temp, nLines, i, i, animIndex, i)
        binIndex = 2
        binIndex(1:nFounders) = 1
        binIndex = binIndex(animIndex)
        temp((nFounders+1):(2*nFounders)) = (/( id, id = 1, nFounders )/)
        animIndex = temp(animIndex)

        ! sorting gamIndex
        call sortix(i, gamIndex, temp) !temp = index for sorting
        gamIndex = gamIndex(temp(1:i))
        j = 1
        if (k .eq. 0) then
           do id = 1, i
              do while(j .ne. gamIndex(id))
                 read(iun,*)
                 j = j + 1
              end do
              read(iun, *) (genome(ichr)%genotypes(animIndex(id), binIndex(id),&
                   iblock) ,iblock = 1, nblock)
              j = j + 1
           end do
        else
           write(STDERR, '(a)') "Sorry"
           write(STDERR, *) "Implementation of this case was removed"
           stop 3
        end if
        close(iun)

        write(filename, '(a,i0,a)') "genome(",iChr,")%positions"
        call alloc1D(genome(iChr)%positions, genome(iChr)%nLoci, &
          trim(filename), "initialiseGenotypes")
        write(filename, '(a,i3.3)') trim(prefixfilename2), ichr
        open(newunit = iun, file = filename, status = 'old')
        read(iun, *)
        do i = 1, genome(ichr)%nloci
           read(iun, *) k , rand
           if ((rand < ZERO) .or. rand > (genome(ichr)%chrL)) then
              write(STDERR, '(a)') "ERROR:"
              write(STDERR, *) " wrong position", rand
              write(STDERR, *) " should be between 0 and ", genome(ichr)%chrL
              stop 2
           end if
           genome(ichr)%positions(i) = rand
        end do
        close(iun)
     end do chr
  end if

  if (verbose) write(STDOUT,'(a)') " end initialising genotypes"

end subroutine InitialiseGenotypes
