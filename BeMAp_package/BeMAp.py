#!/usr/bin/python

import os
import sys
import shutil
import tempfile
import argparse
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from identification import find_AMR, make_each_fasta
from summary import make_summary
from alignment import make_target_centered, align
from mapping import make_dendrogram, make_legend, original_mapping, property_color, property_mapping

parser = argparse.ArgumentParser(description = 'Mapping of the antimicrobial resistance genes in plasmids')

parser.add_argument('-d',
                    '--indir', 
                    required=True, 
                    help = 'directry containing genbank files')
                    
parser.add_argument('-i', 
                    '--infile', 
                    required=True, 
                    help= 'fasta file for identifying the target genes')
                    
parser.add_argument('-o', 
                    '--out', 
                    required=False, 
                    help = 'directory for output')
                    
parser.add_argument('--precise', 
                    required = False,
                    action = 'store_true', 
                    help = 'identify the antimicrobial resistance genes precisely')
                    
parser.add_argument('--ident_target',
                    required = False,
                    type = int,
                    default = 95,
                    help = 'change identity when quering target gene if necessary (default=95%)')
                    
parser.add_argument('--property', 
                    action = 'store_true', 
                    help = 'mapping with Inc group, country, organism and host if necessary') 
                    
parser.add_argument('--num_process', 
                    type = int, 
                    default = 3, 
                    help = 'the number of process (default=3)')

parser.add_argument('--save_all',
                    action = 'store_true',
                    help = 'save all csv files if necessary')
                    
parser.add_argument('--skip_ident',
                    required = False,
                    help = 'Input directory containing csv files after identification of AMR genes and you can skip step for identification')

parser.add_argument('--image', 
                    action = 'store_true', 
                    help = 'output images of mapping if necessary (not work yet)') 

args = parser.parse_args()

logger = logging.getLogger('LogBeMAp')
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
logging.basicConfig(level=logging.INFO, filename='logger.log', format=fmt)

# Types of the antimicrobial resistance genes
gene = pd.DataFrame.from_dict(
 {0: {0: 'target'},
  1: {0: 'aminoglycoside'},
  2: {0: 'beta-lactam'},
  3: {0: 'quinolone'},
  4: {0: 'macrolide'},
  5: {0: 'phenicol'},
  6: {0: 'tetracycline'},
  7: {0: 'trimethoprim'},
  8: {0: 'sulphonamide'},
  9: {0: 'rifampin'},
  10: {0: 'fosfomycin'},
  11: {0: 'polymxycin'},
  12: {0: 'nitrofuran'}},
  orient = 'index')

def genbank_dir(temp, Args):
    if os.path.exists(temp):
        shutil.rmtree(temp)
    if os.path.isdir(Args):
        shutil.copytree(Args, temp)
    else:
        if os.path.isdir(Args):
            shutil.copytree(Args, temp)

def make_fasta_files(temp_fa, temp_gb):
    accession_fa = []
    genbank_fa = []
    os.makedirs(temp_fa, exist_ok = True)
    for i in os.listdir(temp_gb):
        if not i.startswith('.'):
            record = SeqIO.read(temp_gb + i, 'genbank')
            accession_fa.append(record.id)
            genbank_fa.append(i)
            seq = SeqRecord(record.seq, id = record.id, description = record.description)
            SeqIO.write(seq, temp_fa + record.id + '.fasta', 'fasta')
    return accession_fa ,genbank_fa

def multi(x):
    if args.precise:
        result = find_AMR.find_AMR_detail(accession_CDS, gene.loc[x+1,0], x+2, temp1)
    else:
        result = find_AMR.find_AMR(accession_CDS, gene.loc[x+1, 0], x+2, temp1)
    logging.info('-- ' + gene.loc[x+1, 0] + ' resistance gene **** done')
    return [gene.loc[x+1, 0] , result]

if __name__ == '__main__':
    ##### check if outdir exists
    if args.out:
        os.path.isdir(args.out)
    
    ##### make blast database for taget gene #####
    logging.info('start to make blastdb for target gene')
    os.makedirs('database/target', exist_ok = True)
    shutil.copyfile(args.infile, 'database/target/target.fsa')
    blastn_db = 'makeblastdb -in database/target/target.fsa -dbtype nucl'
    os.system(blastn_db)
    
    ##### preparation for analysis #####
    temp_dir1 = tempfile.TemporaryDirectory()
    genbank_dir(temp_dir1.name + '/genbank/', args.indir)
    
    pre_accession ,pre_genbank = make_fasta_files(temp_dir1.name + '/fasta/',temp_dir1.name + '/genbank/')

    accession, genbank = make_summary.find_plasmid(pre_accession, pre_genbank, temp_dir1.name)
    
    ##### make each fasta files #####
    accession_number = []
    os.makedirs(temp_dir1.name + '/each_fasta', exist_ok = True)
    for i in range(len(accession)):
        os.makedirs(temp_dir1.name + '/each_fasta/' + accession[i], exist_ok = True)
        acc_No = make_each_fasta.make_each_fasta(accession[i], 
                                                 temp_dir1.name + '/genbank/' + genbank[i], 
                                                 temp_dir1.name + '/fasta/' + accession[i] + '.fasta',
                                                 temp_dir1.name)
        if acc_No != None:
            accession_number.append(acc_No)
    
    accession_CDS = pd.DataFrame(accession_number, columns = ['accession', 'No_genes'])
    
    if args.skip_ident:
        ##### skip identification #####
        logging.info('skip identification using ' + args.skip_ident)
        AMRs = []
        for i in range(13):       
            AMRs.append(pd.read_csv(args.skip_ident + '/' + str(i) + '.csv',index_col=0))
    else:
        ##### find genes #####
        logging.info('start to find target and AMR genes')
        AMRs = {}
    
        ##### find target gene #####
        AMRs['target'] = find_AMR.find_target(accession_CDS,temp_dir1.name, args.ident_target)
        logging.info('-- target gene **** done')
        
        ##### find AMR genes #####
        temp1 = temp_dir1.name
        p = Pool(args.num_process)
        result = p.map(multi, list(gene.index)[:12])
        p.close()
        p.join()
            
        for each in result:
            AMRs[each[0]] = each[1]
        
        if list(gene[0]) == list(AMRs.keys()):
            AMRs = list(AMRs.values())
            if args.save_all:
                os.makedirs('AMRs', exist_ok = True)
                if args.out:
                    for i in range(len(AMRs)):
                        AMRs[i].to_csv(args.out + 'AMRs/' + str(i) + '.csv')
                else:
                    for i in range(len(AMRs)):
                        AMRs[i].to_csv('AMRs/' + str(i) + '.csv')
        else:
            preAMR = pd.Series(AMRs)
            AMRs = list(preAMR.reindex(index = gene.index).values)
    
    ##### summarize propeties #####
    logging.info('summarize datasets to csv files')
    result_contain_target = make_summary.find_data_containing_target(AMRs[0].astype(int))
    summary = make_summary.summarize_plasmids(pre_accession, pre_genbank, accession, result_contain_target, temp_dir1.name)
    if args.out:
        summary.to_csv(args.out + '/summary.csv')
    else:
        summary.to_csv('summary.csv')
        
    ##### make the target gene centered #####
    logging.info('start alignment')
    target_centered, gene_centered = make_target_centered.make_target_centered(AMRs, temp_dir1.name)
    
    ##### alignment #####
    if align.check_others_filled(target_centered):
        aligned_acc, plasmid_numbering = align.align_target_centered(align.del_filled_lane(target_centered))
    else:
        aligned_acc, plasmid_numbering = align.align_target_centered(target_centered)
    
    if args.save_all:
        if args.out:
            aligned_acc.to_csv(args.out + 'aligned_acc.csv')
            plasmid_numbering.to_csv(args.out + 'plasmid_numbering.csv')
        else:
            aligned_acc.to_csv('aligned_acc.csv')
            plasmid_numbering.to_csv('plasmid_numbering.csv')
            
    ##### make dendrogram #####
    logging.info('make dendrogram and legend')
    temp_dir_figure = tempfile.TemporaryDirectory()
    make_dendrogram.make_dendrogram(plasmid_numbering, temp_dir_figure.name)
    
    ##### make legend #####
    color = pd.DataFrame.from_dict(
    {0: {0: 'none', 1: 'F5F5F5'},
     1: {0: 'target', 1: '000B00'},
     2: {0: 'Aminoglycoside', 1: 'E2421F'},
     3: {0: 'Beta-lactam', 1: '007655'},
     4: {0: 'Quinolone', 1: 'FDD876'},
     5: {0: 'Macrolide', 1: 'E98035'},
     6: {0: 'Phenicol', 1: '6D3C14'},
     7: {0: 'Tetracycline', 1: 'E24F93'},
     8: {0: 'Trimethoprim', 1: '005B98'},
     9: {0: 'Sulfonamide', 1: '008DBD'},
     10: {0: 'Rifampin', 1: '737373'},
     11: {0: 'Fosfomycin', 1: '6846A5'},
     12: {0: 'Polymxycin', 1: 'ABC900'},
     13: {0: 'Nitrofuran', 1: '005E15'}},
     orient = 'index')
    
    color.loc[0,0] = 'Not resistance gene'
    color.loc[3,0] = '$\u03b2$-lactam'
    
    make_legend.make_AMR_legend(color, temp_dir_figure.name)
    
    if args.property:    
        inc_color = property_color.make_inc_color(summary)
        country_color = property_color.make_country_color(summary)
        organism_color = property_color.make_organism_color(summary)
    
        make_legend.make_property_legend('Inc group', inc_color, temp_dir_figure.name)
        make_legend.make_property_legend('country', country_color, temp_dir_figure.name)
        make_legend.make_property_legend('organism', organism_color, temp_dir_figure.name)
    
    temp_dir_figure.cleanup()
    
    ##### make mapping #####
    logging.info('start mapping')
    target_centered_to_int = target_centered.dropna(how='all').fillna(1000).astype(int)
    target_centered_to_int[target_centered_to_int==1000] = np.nan
    new_AMR = pd.DataFrame.reindex(target_centered_to_int, columns=list(aligned_acc.index)).dropna(how='all',axis=0)
    new_centered_gene = pd.DataFrame.reindex(gene_centered, columns=list(aligned_acc.index)).dropna(how='all',axis=0)    
    
    original_mapping.make_original_mapping(new_AMR, new_centered_gene, color, summary, 'figure/legend/AMR.png', 'figure/dendrogram.png', args.out)
    
    if args.save_all:
        if args.out:
            new_AMR.to_csv(args.out + 'new_AMR.csv')
            new_centered_gene.to_csv(args.out + 'new_centered_gene.csv')
        else:
            new_AMR.to_csv('new_AMR.csv')
            new_centered_gene.to_csv('new_centered_gene.csv')

    if args.property:
        color = color.reindex(index=[1,2,3,4,5,6,7,8,9,11,12,13])                                                                       # avoid coloring Not resistance genes
        property_mapping.make_mapping(new_AMR, new_centered_gene, color, inc_color, summary, 
                                      'Inc group', 'figure/legend/AMR.png', 'figure/legend/Inc group.png', 'figure/dendrogram.png', args.out)
        property_mapping.make_mapping(new_AMR, new_centered_gene, color, country_color, summary, 
                                      'country', 'figure/legend/AMR.png', 'figure/legend/country.png', 'figure/dendrogram.png', args.out)
        property_mapping.make_mapping(new_AMR, new_centered_gene, color, organism_color, summary, 
                                      'organism', 'figure/legend/AMR.png', 'figure/legend/organism.png', 'figure/dendrogram.png', args.out)
    
    temp_dir1.cleanup()
    
    logging.info('done')
