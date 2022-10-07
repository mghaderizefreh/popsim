# Pop-sim
This program calculates Gmatrix, simulates breeding values and phenotypes modelled with random regression for a population of sheep given a pedigree structure, correlation structure between breeding values and genepool. A genepool data is available in `/exports/cmvm/eddie/eb/groups/wilson_group/masoud/data/genepool.tar.gz`. Decompress the file and place them in data folder to run the included example. Some of these outputs can be turned off (e.g., Gmatrix). See the example input in test folder.

## Installation
A simple bash script is acompanied with the code. Under eddie load the intel library by

``` shell
module load intel
```

and run the script. The cmake project is currently broken and needs fixing. 

## Usage
Three input files are needed.
1. Input file for all the parameters. See an example in `tests` folder.
2. Pedigree file containing three columns (id, sire, dam). It is **VERY IMPORTANT** that the pedigree is sorted, i.e., each individual identification comes after its parent identification. For those individuals that have unknown parent, 0 must be passed as id of parents. The name of the pedigree file needs to be passed in the input file. There is a simple python script in the folder `scripts` to create a balanced, non-overlapping, full-sib only pedigree structures.
3. Genotype of the base population and SNP positions (one file per chromosome). See an example in `/exports/cmvm/eddie/eb/groups/wilson_group/masoud/data/genepool.tar.gz`. The name of these files must be written in the input file.

To run the program use

``` shell
<path_to_executable> <path_to_input_file> > <path_to_screen_output>
```

If chosen to output the genotype files, they will be written to genotype.txt (or genotype.bin) in the folder from which the above command is executed.

### Run with OMP
The Gmatrix routine is written with OMP loop. If more than one core are available, you can specifiy the number of threads by
``` shell
OMP_NUM_THREADS=<NUMBERS> <path_to_exe> <path_to_input_file> > <path_to_screen_output>
```
Some of small arrays are stored in stack. However, if the number of animals in the pedigree is extremely high stack size must be increased, i.e, prepend the above line with
```shell
KMP_STACKSIZE=<SIZE_in_M_or_G>
```
