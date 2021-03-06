# BeMAp
BeMAp is a tool for mapping and analyzing plasmids carrying the antimicrobial resistance genes on a spreadsheet.

# Features
BeMAp, Bird's-eye MApping of plasmids, is a tool for mappping and analyzing multiple plasmids carrying the antimicrobial resistance genes on a spreadsheet.
Outlining flowchart of BeMAp is below.

<img src="https://user-images.githubusercontent.com/89430651/138423530-dfa53397-c8b3-4dca-aa56-9d4c8f7a6fed.png" width="400">

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
Please download **BeMAp_package**, including all codes to run BeMAp.

If you need test run, you should download **testrun** directory, containing 33 genbank files and fasta files encoding a sequence of blaIMP-6.  
For beginners, testrun directory is recommended to be stored in BeMAp directory.

# Usage
* change directory
```bash
cd BeMAp_package
```

* before perfoming BeMAp, make databases for identification of AMR genes and classification of Inc group
```bash
python makeblastdb.py
```

* testrun (mapping 33 plasmids carrying blaIMP-6) : If you store testrun directory in BeMAp directory
```bash
python BeMAp.py -d testrun/genbank/ -i testrun/testrun.fsa
```

* Now let's start BeMAp!
```bash
python BeMAp.py -d <directory_of_genbank_files> -i <fasta_file_of_target_gene>
```

# Options
```
 -d, --indir [X]     directry containing genbank files 
 -i, --infile [X]    a fasta file for identifying the target gene 
 -o, --out           directory for output 
 --precise           identify the antimicrobial resistance genes precisely 
 --ident_target [X]  change identity when quering target gene if necessary (default=95%) 
 --property          mapping with Inc group, country and organism if necessary 
 --num_process [X]   the number of process (default=3) 
 --save_all          save all csv files if necessary 
 --skip_ident [X]    Input directory containing csv files after identification of AMR genes and you can skip step for identification 
```
Perform BeMAp with precise identification of AMR genes and store all of csv files for analysis
```bash
python BeMAp.py -d directory_of_genbank_files -i fasta_file_of_target_gene --precise --save_all
```
Perform BeMAp using csv files where AMR genes have been identified in BeMAp directory 
```bash
python BeMAp.py -d directory_of_genbank_files -i fasta_file_of_target_gene --skip_ident AMRs/
```

# Databases in BeMAp
BeMAp uses databases for the antimicrobial resistance genes and Inc grouping.
* ResFinder (Final update 4.8.2021)
* PlasmidFinder (Final update 12.7.2021)

# Note
* Displaying the dendric tree depends on PC setting. Please adjust dendric tree to the cells in the spreadsheet.
* Identification of AMR genes precisely takes more times to align multiple plasmids than ambigous identification.
* If you want to download multiple genbank files, use this programm (https://github.com/yusuketsuda/collect_genbank).

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
in submission

# Author
* Yusuke Tsuda
* Bacteriology, Nagoya University Graduate School of Medicine
* tsuda.yusuke(at)med.nagoya-u.ac.jp
