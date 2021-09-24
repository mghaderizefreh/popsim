program simulate
  use constants
  use math
  use global_module
  use rng_module
  use evolution_module
  use blup_module

  implicit none

  logical :: verbose
  logical :: saveSNP, saveQTL
  integer :: nBasePop
  integer :: nanim, nloci, nblock, maxloci, maxblock, ifail, nChr
  character(len=100) :: startFile, filename1, filename2, inputfile
  character(len=100) :: baseNameFreq, gmatrixFile, pedigreeFile, TBVFile
  !integer :: iunlog, iunFarm, iuniter
  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1
  integer, dimension(:), allocatable :: seed
  real(KINDR) :: chrL
  integer, dimension(:), allocatable :: indiv!, male, female
  !  logical, dimension(:), allocatable :: sex !true = male, false = female

  real(KINDR), dimension(:,:), allocatable :: corrs

!  integer, dimension(:,:,:), pointer :: parentGen, offGen

  logical :: randomQTL
  integer :: nfounders
  real(KINDR) :: maf
  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer :: nSNP, nQTL, nComp
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist
  real(KINDR) :: mutationCumP0, chiasmaCumP0, mutationRate
  integer :: maxchiasma, maxmutations
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP, mutationCumP
  integer, allocatable, dimension(:) :: chr_nlocibefore
  integer, allocatable, dimension(:) :: totalChiasma, totalMutation
  real(KINDR), allocatable, dimension(:,:) :: tbv
  real(KINDR), allocatable, dimension(:) :: Amat
  integer, allocatable, dimension(:,:) :: pedigree
  integer :: ivar, iscaled!, maxid
  ! counters
  integer :: i, j, k, iChr, id, igam, iparent
  character(len=30) :: status, estatus

  ! ==============================================
  ! initialising seed
  ! ==============================================
  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(STDERR, 100) "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  ! ==============================================
  ! reading input file
  ! ==============================================
  status = "old"
  call askFileName(inputfile, "input file: ", status, estatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, 101) "error in openning input file ", trim(inputfile)
     stop 2
  end if

  call readInput(inputfile, verbose, nchr, filename1, filename2, nBasePop,&
       chrL, mutationRate, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, &
       corrs, nanim, pedigreefile, gMatrixfile, TBVfile, saveQTL, saveSNP)

  ! ==============================================
  ! allocations
  ! ==============================================
  call alloc1I(indiv, nanim, "indiv", "main")
  indiv= (/( i, i = 1, nanim )/)
  allocate(genome1(nChr))
  call alloc2I(SNPlist, nChr, nSNP, "SNPlist", "main")
  call alloc2I(QTLlist%indices, nChr, nQTL, "QTLlist%indices", "main")
  call alloc3D(QTLlist%values, nChr, nQTL, nComp, "QTLlist%values", "main")
  call alloc2D(tbv, nanim, ncomp, "tbv", "main")
  i = nAnim * (nAnim + 1) / 2 
  call alloc1D(AMat, i, "AMat", "main")
  call alloc1I(chr_nlocibefore, nchr, "chr_nlociBefore", "main")

  ! ==============================================
  ! reading pedigree file
  ! ==============================================
  call alloc2I(pedigree, nanim, 3, "pedigree", "main")
  call readPedigree(pedigreeFile, nanim, pedigree, nfounders)

  ! ==============================================
  ! Preparing genome
  ! ==============================================
  if (nBasePop < nfounders) then 
     write(STDERR, '(a)') "Sorry"
     write(STDERR, *) "There are too many founders in the pedigree file"
     write(STDERR, *) "The program cannot continue. Data files must be expanded"
     stop 3
  end if
  ! ----------------------------------------------
  ! calculating the cumulative distribution for number of recombination events
  ! ----------------------------------------------
  i = 200 ! upto 200 recombinations
  if (verbose) write(STDOUT, 100) &
       "Setting total number of mutations and recombinations"
  call GetMutRecArray(verbose, i, chrL, mutationRate, nLoci, &
       chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma, &
       mutationCumP, mutationCumP0, totalMutation, maxmutations)
  if (verbose) write(STDOUT, 100) "Mutation and recombinations set"

  i = 3 ! genstart = 3
  j = 1 ! istore = 1
  call InitialiseGenotypes(verbose, nchr, nAnim, i, nloci, nblock, j, &
       nfounders, genome1, maxloci, maxblock, ifail, chrL, trim(filename1),&
       trim(filename2))
  ! now adding other animals
  i = 1 ! istore = 1
  k = 1 ! samepos = 1
  do ichr = 1, nChr
     do id = nFounders + 1, nAnim
        do igam = 1, 2
           iparent = pedigree(id, igam + 1)
           call SampleGamete(iparent, id, igam, genome1(iChr)%genotypes,&
                genome1(iChr)%nLoci, genome1(iChr)%nBlock, maxChiasma, &
                chiasmaCumP0, ChiasmaCumP, i, &
                offGenInp = genome1(iChr)%genotypes,&
                positions = genome1(iChr)%positions, samePos = k)
           call sampleMutation(id, igam, genome1(iChr)%genotypes, &
                genome1(iChr)%nLoci, maxMutations, mutationCumP0, &
                mutationCumP, i)
        end do
     end do
  end do

  if (verbose) write(STDOUT, 100) "Genome ready"

  ! ================================================================
  ! Getting QTL and SNP list
  ! ================================================================
  if (verbose) write(STDOUT, 100) "Getting SNPs and their values"
  call getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, randomQTL, genome1,&
       QTLlist, SNPlist, corrs, baseNameFreq, maf)
  if (verbose) write(STDOUT, 100) "QTL and SNP list simulated"
  ! saving QTL list (if instructed)
  if (saveQTL) then
     open(1, file = "QTLlist.txt")
     do ichr = 1, nchr
        do i = 1, nQTL
           write(1, *) QTLlist%values(ichr, i, 1:ncomp)
        end do
     end do
     close(1)
  end if
  ! saving SNP list (if instructed)
  if (saveSNP) then
     open(1, file = "SNPlist.txt")
     do ichr = 1, nchr
        do i = 1, nSNP
           write(1, '(2i8)') ichr, SNPlist(ichr, i)
        end do
     end do
     close(1)
  end if
  
  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  if (verbose) write(STDOUT, 206) "simulating BV for generation", 0
  call SimulateTBV(nAnim, nChr, nComp, indiv, genome1, chr_nlocibefore,&
       QTLlist, TBV, verbose)
  if (verbose) write(STDOUT, 100) "breeding values simulated"
  
  ! saving TBV (one of outputs)
  open(1,file = trim(tbvFile))
  do i = 1, nanim
     write(1, *) (tbv(i,j), j = 1, ncomp)
  end do
  close(1)

  ! ================================================================
  ! making G matrix
  ! ================================================================
  iscaled = 1 !(0:no, 1:yes)
  ivar  = 1 !(0:sample, 1:2pq, 2:2p'q')
  j = 0 ! imiss (0:mean, 1:ignore)
  i = 1 ! addDom (1:additive, 2:dominance)
  call getGmatrix(nanim, nChr, nSNP, indiv, genome1, SNPlist,&
       chr_nlocibefore, iscaled, ivar, j, i, Amat, verbose)
  !writing AMAT to disk (one of outputs)
  i = maxval(indiv)
  k = 1
  do while (i >= 10)
     i = i / 10
     k = k + 1
  end do
  open(1, FILE = trim(GmatrixFile))
  i = nanim * ( nanim + 1 ) / 2
  write(formato,'(a2,i1,a5,i1,a22)')'(i',k,',1x,i',k,',1x,g24.15,i9,g24.15)'
  write(6,*) " formato= ", trim(formato)
  k = 0
  do i = 1, nanim
     do j = 1, i
        k = k + 1
        write(1, formato) indiv(i), indiv(j), amat(k)
     end do
  end do
  close(1)

  ! ================================================================
  ! setting a new seed for next run
  ! ================================================================
  call ifinal(seed, startfile)

  ! ================================================================
  ! formats
  ! ================================================================
100 format(a)
101 format(2a)
206 format(a, i2)

end program simulate
