# BeMAp
BeMAp is a tool for mapping and analyzing plasmids carrying antimicrobial resistance genes on a spreadsheet.

# Features
BeMAp, Bird's-eye MApping of plasmids, is a tool for mappping and analyzing multiple plasmids carrying antimicrobial resistance genes on a spreadsheet.
Outlining flowchart of BeMAp is below.

![Fig1-11](https://user-images.githubusercontent.com/89430651/138423530-dfa53397-c8b3-4dca-aa56-9d4c8f7a6fed.png)

If you need more detail, please check our article or send e-mail.


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
Please download **BeMAp.package**. This includes all codes to run BeMAp.


If you need test run, you should download **testrun** directory.

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
|* -d, --indir|directry containing genbank files|
|* -i, --infile|a fasta file for identifying the target gene|
* -o, --out
directory for output
* --precise
identify the antimicrobial resistance genes precisely
* --ident_target
change identity when quering target gene if necessary (default=95%)
* --property
mapping with Inc group, country and organism if necessary
* --num_process
the number of process (default=3)
* --save_all
save all csv files if necessary
* --skip_ident
Input directory containing csv files after identification of AMR genes and you can skip step for identification

# Databases in BeMAp
BeMAp uses databases for the antimicrobial resistance genes and Inc grouping.
* ResFinder (Final update)
* PlasmidFinder (Final update)

# Note
* Displaying the dendric tree depends on PC setting. Please adjust dendric tree to the cells in the spreadsheet.

# Future Plan
BeMAp will offer options below in the future.

* Provide images of mapping not a spreadsheet for publication or presentation
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
