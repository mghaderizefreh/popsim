subroutine readInput(inputfile, verbose, nchr, genepoolfile, geneposfile, &
     nBasePop, chrL, mu, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, &
     corrs, nanim, pedigreefile, gMatFile, tbvFile, saveQTL, saveSNP)
  use constants
  implicit none
  character(len=*), intent(in) :: inputfile ! list of all inputs
  logical, intent(out) :: verbose
  !!!!!!!!!!!!!!!!!!!! genomic !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nchr
  character(len=100), intent(out) :: genepoolfile
  character(len=100), intent(out) :: geneposfile
  integer, intent(out) :: nbasePop
  real(KINDR), intent(out) :: chrL
  real(KINDR), intent(out) :: mu
  integer, intent(out) :: nQTL
  integer, intent(out) :: nSNP
  logical, intent(out) :: randomQTL
  real(KINDR), intent(out) :: MAF
  character(len=100), intent(out) :: baseNameFreq
  !!!!!!!!!!!!!!!!!!!! genetic !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nComp
  real(KINDR), allocatable, dimension(:,:), intent(out) :: corrs
  !!!!!!!!!!!!!!!!!!!! population !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nanim
  !!!!!!!!!!!!!!!!!!!! output !!!!!!!!!!!!!!!!!!!!
  character(len=100), intent(out) :: pedigreeFile
  character(len=100), intent(out) :: gMatFile
  character(len=100), intent(out) :: tbvFile
  !!!!!!!!!!!!!!!!!!!! debugging files !!!!!!!!!!!!!!!!!!!!
  logical, intent(out) :: saveQTL
  logical, intent(out) :: saveSNP
  !!!!!!!!!!!!!!!!!!!! dummy variables !!!!!!!!!!!!!!!!!!!!
  character(len = 200) :: line, formato
  integer :: iinput, stat, iun, lno, j, i
  real(KINDR) :: rinput
  logical :: bool

  lno = 0
  open(newUnit = iun, file = trim(inputfile))

33 format(a34,": ", l1)
34 format(a34,": ", i0)
35 format(a34,": ", g0.15)
36 format(a34,": ", a)

  ! verbosity
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for verbose", lno)
  verbose = iinput .eq. 1
  write(STDOUT, 33) "verbose?", verbose

  !!!!!!!!!!!!!!!!!!!! genomic !!!!!!!!!!!!!!!!!!!!
  ! number of chrosomomes
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for nchr", lno)
  call assert(iinput.gt.0, "nchr must be greater than 0", lno)
  nchr = iinput
  write(STDOUT, 34) "number of chromosomes", nchr

  ! base-name for genepool filename
  call nextInput(iun, line, lno)
  genepoolfile = trim(line)
  do i = 1, nchr
     write(line, '(a,i3.3)') trim(genepoolfile),i
     inquire(file=line, exist = bool)
     write(formato, '(a12,a)') "cannot find ", trim(line)
     call assert(bool, formato, lno)
  end do
  write(STDOUT, 36) "file for genepool", genepoolfile

  ! n of animals in genepool (using `line`  [don't remove the previousblock])
  open(1, file = trim(line))
  read(1, *) nBasePop
  write(STDOUT, 34) "number of animals in genepool", nBasePop
  close(1)

  ! base-name for position filename
  call nextInput(iun, line, lno)
  geneposfile = trim(line)
  do i = 1, nchr
     write(line, '(a,i3.3)') trim(geneposfile), i
     inquire(file=line, exist = bool)
     call assert(bool, &
          "error in finding a file for genepos with as many as chromosomes", lno)
  end do
  write(STDOUT, 36) "file for genepositions", geneposfile

  ! chrL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read input for chrL", lno)
  call assert(rinput.gt.ZERO, "Chromosome length must be > 0.0", lno)
  chrL = rinput
  write(STDOUT, 35) "chromosome length", chrL

  ! mu (mutation rate)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read input for mu", lno)
  call assert(rinput.gt.ZERO, "mutation rate must be > 0.0", lno)
  mu = rinput
  write(STDOUT, 35) "mutation rate", mu

  ! nQTL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of QTLs", lno)
  call assert(iinput.gt.0, "nQTL must be > 0", lno)
  nQTL = iinput
  write(STDOUT, 34) "number of QTL/chromosome", nQTL

  ! nSNP
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of SNPs", lno)
  call assert(iinput.gt.0, "number of SNP cannot be < 0", lno)
  nSNP = iinput
  write(STDOUT, 34) "number of SNP/chromosome", nSNP

  ! randomQTL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for randomQTL", lno)
  randomQTL = iinput .eq. 1
  write(STDOUT, 33) "Are QTLs random?", randomQTL

  ! MAF
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read MAF for finding QTLs", lno)
  if (.not.randomQTL) then
     call assert((rinput.ge.ZERO).and.(rinput.lt.HALF), "maf must be &
          &between 0.0 and 0.5 (incl,excl., respecitvely)",lno)
     maf = rinput
     write(STDOUT, 35) "MAF cutoff for QTL", MAF
  else
     maf = rinput
     write(STDOUT, 35) "MAF cutoff for QTL (is ignored)", MAF
  end if

  ! baseNameFreq
  call nextInput(iun, line, lno)
  baseNameFreq = trim(line)
  if (.not.randomQTL) then
     do i = 1, nchr
        write(line, '(a,i3.3)') trim(baseNameFreq),i
        inquire(file=line, exist = bool)
        call assert(bool, "error in finding a file for frequenecies with&
             &as many as chromosomes", lno)
     end do
     write(STDOUT, 36) "file for frequnecy", baseNameFreq
  else
     write(STDOUT, 36) "file for frequency (is ignored)", baseNameFreq
  end if

  !!!!!!!!!!!!!!!!!!!! genetic !!!!!!!!!!!!!!!!!!!!
  ! nComp
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of components", lno)
  ncomp = iinput
  write(STDOUT, 34) "number of components", nComp

  ! correlations
  call alloc2D(corrs, ncomp, ncomp, "corrs", "readInput")
  if (ncomp > 1) then
     do i = 2, ncomp
        do j = 1, (i - 1)
           call nextInput(iun, line, lno)
           read(line, *, iostat = stat) rinput
           call assert(stat.eq.0, "failed to read a correlation", lno)
           call assert((rinput.ge.-1._KINDR).and.(rinput.le.ONE),&
                "correlation must be in [-1,+1]", lno)
           corrs(i,j) = rinput
           corrs(j,i) = rinput
           write(formato, '(a5,2(i1,a1))') "corr(", i, ",", j, ")"
           write(STDOUT, 35) trim(formato), corrs(i,j)
        end do
     end do
  end if
  do i = 1, ncomp
     corrs(i,i) = ONE
  end do
  !!!!!!!!!!!!!!!!!!!! population !!!!!!!!!!!!!!!!!!!!
  ! nanim
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for nanim", lno)
  nanim = iinput
  write(STDOUT, 34) "number of animals", nanim  

  ! pedigreeFile
  call nextInput(iun, line, lno)
  pedigreeFile = trim(line)
  inquire(file=pedigreeFile, exist = bool)
  call assert(bool, "could not find pedigree file", lno)
  write(STDOUT, 36) "pedigree filename", trim(pedigreeFile)
  !!!!!!!!!!!!!!!!!!!! output !!!!!!!!!!!!!!!!!!!!
  ! gMatFile
  call nextInput(iun, line, lno)
  gMatFile = trim(line)
  write(STDOUT, 36) "Gmatrix filename", trim(gMatFile)

  ! TBVfile
  call nextInput(iun, line, lno)
  TBVfile = trim(line)
  write(STDOUT, 36) "TBV filename", trim(TBVfile)
  !!!!!!!!!!!!!!!!! debuggin files !!!!!!!!!!!!!!!!!
  ! SaveQTL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for SaveQTL", lno)
  saveQTL = iinput .eq. 1
  write(STDOUT, 33) "saveQTL?", saveQTL

  ! SaveSNP
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for SaveSNP", lno)
  saveSNP = iinput .eq. 1
  write(STDOUT, 33) "saveSNP?", saveSNP
  !!!!!!!!!!!!!!!!!!!! FINISHED !!!!!!!!!!!!!!!!!!!!
  write(STDOUT, '(a,i4,a)') "all inputs read in", lno, " lines"
  close(iun)
end subroutine readInput
subroutine nextInput(iun, output, lno)
  use constants
  implicit none
  character(len = 200), intent(out) :: output
  integer, intent(in) :: iun
  integer, intent(inout) :: lno
  character(len = 200) :: line
  do
     read(iun, '(A)', end = 100) line
     lno = lno + 1
     if(line(1:1) .eq. '!') cycle
     if(trim(line) .eq. '') cycle
     output = trim(line)
     return
100  write(STDERR, *) "ERROR: unexpected end of input file" 
     write(STDERR, *) lno, "lines were read"
     stop 2
  end do
end subroutine nextInput
subroutine assert(iostat, msg, no)
  use constants
  implicit none
  character(len = *), intent(in) :: msg
  logical, intent(in) :: iostat
  integer, intent(in) :: no
  if (iostat) return
  write(STDERR, '(a21,i3)') "ERROR in line ", no
  write(STDERR, *) msg
  stop 2
end subroutine assert
