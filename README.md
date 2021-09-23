# BeMAp
BeMAp is a platform for mapping and analyzing plasmids carrying antimicrobial resistance genes on a spreadsheet.

# Features


# Requirements 
* NCBI blast
* Biopython

# Tested environments
* macOS 11.6
* Ubuntu 20.04
* Python 3.7.10, 3.8.5
* makeblastdb 2.5.0, 2.9.0
* biopython 1.77, 1.79

# Installation


# Usage
```bash
cd BeMAp_package
```

* before perfoming BeMAp, make database for AMR genes and Inc typing
```bash
python makeblastdb.py
```

* Now let's start BeMAp!
```bash
python BeMAp.py -d directory_of_genbank_files -i fasta_file_of_target_gene
```

* testrun
```bash
cd BeMAp.package
python BeMAp.py -d sample -i sample.fsa
```

# Options
* -d, --indir 

* -i, --infile

* -o, --out

* -num_process

* -identity



# Note

# License

# Reference


# Author
* Yusuke Tsuda
* Bacteriology, Nagoya University Graduated School of Medicine
* tsuda.yusuke@med.nagoya-u.ac.jp
