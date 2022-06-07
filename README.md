# AssemblyDownGrader

Creates draft-like assemblies by chopping complete reference genomes and possibly introducing mutations. The goal is to obtain draft-like genome sequences without the computationally intensive steps of read simulation and genome assembly. Useful to benchmark bioinformatic tools designed to analyse draft genomes.

To install this tool you need to clone this repository to your local machine and install dependencies using the `yaml` file supplied. 

Dependencies are:

- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [bioawk](https://github.com/lh3/bioawk)
- [Mutation-Simulator](https://github.com/mkpython3/Mutation-Simulator) with python v3.10 installed (need pip to install)

The script is developed in `bash` and can be run from the command line using any GNU/Linux distribution by typing `assemblydowngrader.sh` after adding the location of the script to $PATH variable. This tool can be parametrized from the command line. Each run gets its own unique ID present in the output file names, which ensures that previous data do not get overwritten.

The options are the following:

````
        -i|--infile
        Specifies the input high-quality genome sequence. The only mandatory option.

        -minp|--minpoints
        The minimum number of breakpoints to chop up the input sequences. [default = 1]

        -maxp|--maxpoints
        The minimum number of breakpoints to chop up the input sequences. Within the range specified by -minp and -maxp a random number of breakpoints will be chosen for each contig. [default = 5]
        
        -r|--reps
        Number of independent replications. This many alternative draft-like assemblies will be produced with a unique error profile and random gaps. [default = 3]
        
        -ming|--mingaps
        Minimum length of sequencing gaps introduced. [default = 1]
        
        -maxg|--maxgaps
        Maximum length of sequencing gaps introduced. A random length within the range specified by -ming and -maxg will be assigned to each sequencing gap. [default = 20]

        -minlen|--minimum-length
        Minimum lenght of contigs to simulate. Does not influence the position of random breaks but filters the resulting contigs. [default = 1]

        -out|--out-prefix
        Prefix of the output files. [default = draft]

        -sn|--snp-rate
        The rate at which SNPs should be simulated. [default = 0.0002]

        -in|--insertion-rate
        The rate at which insertions should be simulated. [default = 0.0001]

        -de|--deletion-rate
        The rate at which deletions should be simulated. [default = 0.0001]
     
        -iv|--inversion-rate
        The rate at which inversions should be simulated. [default = 0.00005]

        --inmax|--insert-max-length
        The maximum length of insert mutations. Minimum length is 1. [default = 10]	
        
        --ivmax|--inversion-max-length
        The maximum length of inversions. Minimum length is 1. [default = 100]
        
        --demax|--deletion-max-length
        The maximum length of deletions. Minimum length is 1. [default = 5]

        --du|--duplication-rate
        The rate at which (segmental) duplications should be simulated. [default = 0 meaning no duplications]
        
        --tl|--translocation-rate
        The rate at which translocations should be simulated. [default = 0 meaning no translocations]

        --dumin|--duplication-min-length
        The minimum length of duplications. [default = 1]
 
        --dumax|--duplication-max-length
        The maximum length of duplications. [default = 2]

        --tlmin|--translocation-min-length
        The minimum length of translocations. [default = 1]

        --tlmax|--translocation-max-length
        The maximum length of translocations. [default = 2]
        

**for a more detailed explanation of mutation parameters, please see https://github.com/mkpython3/Mutation-Simulator
````

## **This is a work in progress. Use it only at your own risk and doublecheck your data. Feedbacks are very welcome.**
