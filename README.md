# Pop-sim
This program calculates Gmatrix and simulates breeding values for a population of sheep given a pedigree structure, correlation structure between breeding values and genepool. A genepool data is available in `/exports/cmvm/eddie/eb/groups/wilson_group/masoud/data/genepool.tar.gz`. Untar the file and place them in data folder to run the included example.

## Installation
A simple bash script is acompanied with the code. Under eddie load the intel library by

``` shell
module load intel
```

A cmake project is also available in case you need to debug something. Note you need to load a decent version of gcc and cmake under eddie for it to work.

It is also available to compile with gfortran. 

## Usage
Three input files are needed.
1. Input file for all the parameters. See an example in `tests` folder.
2. Pedigree file containing three columns (id, sire, dam). It is **VERY IMPORTANT** that the pedigree is sorted, i.e., each individual identification comes after its parent identification. For those individuals that have unknown parent, 0 must be passed as id of parents. The name of the pedigree file needs to be passed in the input file.
3. Genotype of the base population and SNP positions (one file per chromosome). See an example in `/exports/cmvm/eddie/eb/groups/wilson_group/masoud/data/genepool.tar.gz`. The name of these files must be written in the input file.

To run the program use

``` shell
echo <name_of_the_input_file> | <path_to_executable> > <path_to_screen_output>
```


