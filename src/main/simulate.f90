program simulate
  use constants
  use math
  use global_module
  use rng_module
  use evolution_module
  use blup_module

  implicit none

  logical :: verbose
  logical :: saveSNP, saveQTL, saveGenotypeBIN, saveGenotypeTXT, saveGRM
  logical :: saveFullGenome
  integer :: nBasePop
  integer :: nanim, nloci, nblock, maxloci, maxblock, ifail, nChr
  character(len=100) :: startFile, filename1, filename2, inputfile, phenFile
  character(len=100) :: baseNameFreq, gmatrixFile, pedigreeFile, TBVFile, genotypeFileBIN, genotypeFileTXT
  !integer :: iunlog, iunFarm, iuniter
  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1
  integer, dimension(:), allocatable :: seed
  real(KINDR) :: chrL
  integer, dimension(:), allocatable :: indiv!, male, female
  type(variances) :: vars
  real(KINDR), dimension(:), allocatable :: means, cv, h2, phen
  real(KINDR) :: xmin, xmax, val1
  integer :: nfarm
  real(KINDR), dimension(:,:), allocatable :: locations
  real(KINDR) :: farmRange
  integer :: allocation, nlox
  real(KINDR), dimension(:,:), allocatable :: corrs
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
  integer, allocatable, dimension(:) :: chr_nlocibefore, farmInd
  integer, allocatable, dimension(:) :: totalChiasma, totalMutation
  real(KINDR), allocatable, dimension(:,:) :: tbv
  real(KINDR), allocatable, dimension(:) :: Amat
  integer, allocatable, dimension(:,:) :: pedigree
  integer :: ivar, iscaled!, maxid
  ! counters
  integer :: i, j, k, iChr, id, igam, iparent
  character(len=30) :: status, estatus
  external :: setUserParameters, simulatePhenotype
  ! ==============================================
  ! initialising seed
  ! ==============================================
  startfile = "seed.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(STDERR, 100) "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  ! ==============================================
  ! reading input file
  ! ==============================================
  estatus = "old"
  i = command_argument_count()
  if (i == 0) then
   call askFileName(inputfile, "input file: ", status, estatus)
   if (status(1:1) .eq. "x") then
     write(STDERR, 101) "error in openning input file ", trim(inputfile)
     stop 2
   end if
  elseif (i >= 1) then
   call get_command_argument(1, inputfile)
  end if

  call readInput(inputfile, verbose, nchr, filename1, filename2, nBasePop,&
       chrL, mutationRate, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, &
       corrs, nanim, pedigreefile, TBVFile, phenFile, saveGRM, gMatrixfile, &
       saveGenotypeBIN, saveGenotypeTXT, genotypeFileBIN, genotypeFileTXT, &
       vars, means, cv, h2,xmin, xmax, nlox, nfarm, farmRange, allocation,&
       saveQTL, saveSNP, saveFullGenome)

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
  i = nanim * nlox
  call alloc2D(locations, nanim, nlox, "x", "main")
  call alloc1I(farmInd, i, "farmInd", "simulatePhenotype")
  call alloc1D(phen, i, "phen", "main")
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
  ! initialising the genotypes (for founders)
  ! ----------------------------------------------
  i = 3 ! genstart = 3
  j = 1 ! istore = 1
  call InitialiseGenotypes(verbose, nchr, nAnim, i, nloci, nblock, j, &
       nfounders, genome1, maxloci, maxblock, ifail, chrL, trim(filename1),&
       trim(filename2))
  ! ----------------------------------------------
  ! calculating the cumulative distribution for number of recombination events
  ! ----------------------------------------------
  i = 200 ! upto 200 recombinations
  if (verbose) write(STDOUT, 100) &
       "Setting total number of mutations and recombinations"
  ! This is needed to calculate "an average" number of mutations per chromosome
  ! even though, the mutation is totally ignored here.
  ! in a better scenario, one should calculate total number of mutations
  ! independently
  nloci = 0
  do iChr = 1, nChr
     nloci = nloci + genome1(iChr)%nLoci
  end do
  nLoci = nLoci / nChr
  call GetMutRecArray(verbose, i, chrL, mutationRate, nLoci, &
       chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma, &
       mutationCumP, mutationCumP0, totalMutation, maxmutations)
  if (verbose) write(STDOUT, 100) "recombinations set (mutations ignored)"

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
     write(formato, '(a1,i2,a17)') "(", ncomp-1, "(f15.7,1x),f15.7)"
     if (verbose) write(STDOUT, *) "formato: ", trim(formato)
     open(1, file = "QTLlist.txt")
     do ichr = 1, nchr
        do i = 1, nQTL
           write(1, formato) QTLlist%values(ichr, i, 1:ncomp)
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

  ! saving whole genome for both haplotypes and all markers (if instructed)
  if (saveFullGenome) then
    do ichr = 1, nchr
      write(formato, '(a,i0.3)') "genome.ch", ichr
        open(1, file = trim(formato))
        write(1, '(2(i0,1x),i0)') nanim, genome1(ichr)%nloci, genome1(ichr)%nblock
        write(formato, '(a, i0, a)') "(", genome1(ichr)%nblock-1, "(i0,1x),i0)"
        do i = 1, nanim
          write(1, fmt=formato) &
              (genome1(ichr)%genotypes(i, 1, j), j = 1, genome1(ichr)%nblock)
          write(1, fmt=formato)&
              (genome1(ichr)%genotypes(i, 2, j), j = 1, genome1(ichr)%nblock)
        end do
        close(1)
    end do
  end if
  
  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  if (verbose) write(STDOUT, 206) "simulating BV for generation", 0
  call SimulateTBV(nAnim, nChr, nComp, indiv, genome1, chr_nlocibefore,&
       QTLlist, TBV, verbose)
  if (verbose) write(STDOUT, 100) "breeding values simulated"
  ! scale 
  do i = 1, nComp
   val1 = sqrt(vars%A(i) / variance(TBV(1:nanim,i), nanim))
   tbv(1:nanim,i) = tbv(1:nanim,i) * val1
   val1 = sum(tbv(1:nanim,i)) / nanim
   tbv(1:nanim,i) = tbv(1:nanim,i) - val1
  end do
  
  ! saving TBV (one of outputs)
  open(1,file = trim(tbvFile))
  write(formato, '(a1,i2,a17)') "(", ncomp-1, "(f15.7,1x),f15.7)"
  if (verbose) write(STDOUT, *) "formato: ", trim(formato)
  do i = 1, nanim
     write(1, formato) (tbv(i,j), j = 1, ncomp)
  end do
  close(1)
  ! ================================================================
  ! making G matrix
  ! ================================================================
  if (saveGRM) then
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
  end if
  ! ================================================================
  ! writing genotype file
  ! ================================================================
  if (saveGenotypeBIN .or. saveGenotypeTXT) then
   call makeGenotype(nanim, nchr, nsnp, indiv, genome1,&
   SNPlist, chr_nlocibefore, genotypeFileBIN, genotypeFileTXT,&
   saveGenotypeBIN, saveGenotypeTXT)
  end if
  ! ================================================================
  ! simulating and writing phenotype
  ! ================================================================  
  i = 2 ! ncomp
  call setUserParameters(xmin, xmax, nfarm, farmRange, allocation)
  call SimulatePhenotype(verbose, nanim, nComp, tbv, means, vars, nlox,&
   phen, locations, farmInd, pedigree)
  open(1, file = trim(phenFile))
  k = 0
  do i = 1, nanim
   do j = 1, nlox
    k = k + 1
    write(1, '(i0,1x,i0,1x, g0.14, 1x, g0.15)') i, farmInd(k), locations(i, j) , phen(k)
   end do
  end do

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
