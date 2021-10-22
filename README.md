# BeMAp
BeMAp is a tool for mapping and analyzing plasmids carrying antimicrobial resistance genes on a spreadsheet.

# Features
BeMAp, Bird's-eye MApping of plasmids, is a tool for mapppin and analyzing multiple plasmids carrying antimicrobial resistance genes on a spreadsheet.
Outlining flowchart of BeMAp is below.




# Requirements 
* NCBI blast
* Biopython

# Tested environments
* macOS 11.6
* Ubuntu 20.04
* Python 3.7.10 and 3.8.5
* makeblastdb 2.5.0 and 2.9.0
* biopython 1.77 and 1.79
* Microsoft Excel for Mac 16.53
* Libre office 6.4.7.2

# Installation


# Usage
* change directory
```bash
cd BeMAp_package
```

* before perfoming BeMAp, make databases for AMR genes and Inc typing
```bash
python makeblastdb.py
```

* testrun (mapping of 20 plasmids carrying blaIMP-6)
```bash
python BeMAp.py -d testrun -i testrun.fsa
```

* Now let's start BeMAp!
```bash
python BeMAp.py -d directory_of_genbank_files -i fasta_file_of_target_gene
```

# Options
* -d, --indir 

* -i, --infile

* -o, --out

* -num_process

* -identity



# Future Plan
BeMAp will offer options below in the future.
* Precise identification of the antimicrobial resistance genes
* Precise identification of the other genes, such as virulence genes or transposases
* New alignments algorithm by weighting according to features of plasmids or organizations of genes
* Analyse including chromosomes
Please wait, thank you.

# License
MIT License

# Reference


# Author
* Yusuke Tsuda
* Bacteriology, Nagoya University Graduated School of Medicine
* tsuda.yusuke@med.nagoya-u.ac.jp
