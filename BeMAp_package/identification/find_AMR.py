#!/usr/bin/python

from Bio import SeqIO
from Bio.Blast import Applications
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
import io
import re

def find_AMR(accession_name_No_genes, reagent, number, temp):
    accession = accession_name_No_genes['accession']
    No_genes = accession_name_No_genes['No_genes']

    # find AMR genes. cells in no gene is to be 0.01 
    detect_AMR = pd.DataFrame(np.array([[0.01]*len(accession) for i in range(max(No_genes))]),columns = accession)
    for i in accession_name_No_genes.index:
        for k in range(accession_name_No_genes.loc[i, 'No_genes']):
            cline = Applications.NcbiblastnCommandline(query = temp + '/each_fasta/' + accession_name_No_genes.loc[i, 'accession'] +'/' + str(k) + '.fasta', db='database/ResFinder_db/' + reagent + '.fsa', outfmt=5)()[0]
            blast_result=NCBIXML.read(io.StringIO(cline))
            description = blast_result.descriptions
            if description != []:
                for l in  description:
                    if l.e == 0:
                        detect_AMR.iloc[k,i] = number
                        break
                    else:
                        detect_AMR.iloc[k,i] = 0
            else:
                detect_AMR.iloc[k,i] = 0

    return detect_AMR


def find_target_old(accession_name_No_genes, temp, identity=100):
    accession = accession_name_No_genes['accession']
    No_genes = accession_name_No_genes['No_genes']

    # find AMR genes. cells in no gene is to be 0.01
    detect_target = pd.DataFrame(np.array([[0.01]*len(accession) for i in range(max(No_genes))]),columns = accession)
    for i in accession_name_No_genes.index:
        for k in range(accession_name_No_genes.loc[i, 'No_genes']):
            cline = Applications.NcbiblastnCommandline(query = temp + '/each_fasta/' + accession_name_No_genes.loc[i, 'accession'] +'/' + str(k) + '.fasta', db='database/target/target.fsa', perc_identity=identity, outfmt=5)()[0]
            blast_result=NCBIXML.read(io.StringIO(cline))
            description = blast_result.descriptions
            if description != []:
                for l in  description:
                    if l.e == 0:
                        detect_target.iloc[k,i] = 1
                        break
                    else:
                        detect_target.iloc[k,i] = 0
            else:
                detect_target.iloc[k,i] = 0

    return detect_target

def find_target(accession_name_No_genes, temp, identity=95):
    '''
    for identification of the antimicrobail resistance genes roughly
    '''
    accession = accession_name_No_genes['accession']
    No_genes = accession_name_No_genes['No_genes']

    # find AMR genes. cells in no gene is to be 0.01
    detect_target = pd.DataFrame(np.array([[0.01]*len(accession) for i in range(max(No_genes))]),columns = accession)
    for i in accession_name_No_genes.index:
        for k in range(accession_name_No_genes.loc[i, 'No_genes']):
            cline = Applications.NcbiblastnCommandline(query = temp + '/each_fasta/' + accession_name_No_genes.loc[i, 'accession'] +'/' + str(k) + '.fasta', db='database/target/target.fsa', perc_identity=identity, outfmt=5)()[0]
            blast_result=NCBIXML.read(io.StringIO(cline))
            description = blast_result.descriptions
            if description != []:
                for l in  description:
                    if l.e == 0:
                        detect_target.iloc[k,i] = 1
                        break
                    else:
                        detect_target.iloc[k,i] = 0
            else:
                detect_target.iloc[k,i] = 0

    return detect_target

def revise_target(pd, temp, identity=100):
    pd_revise = pd.copy()

    num_gene = re.compile(r'\|\|(.*)?')

    for i in pd_revise.index:
        for j in pd_revise.columns:
            if int(pd_revise.loc[i,j]) == 1:
                cline = Applications.NcbiblastnCommandline(query = temp + '/each_fasta/' + str(j) +'/' + str(i) + '.fasta', 
                                                           db='database/ResFinder_db_detail/beta-lactam_detail.fsa', 
                                                           perc_identity=identity, 
                                                           outfmt=5)()[0]
                blast_result=NCBIXML.read(io.StringIO(cline))
                description = blast_result.descriptions
                if description != []:
                    for l in  description:
                        if l.e == 0:
                            pd_revise.loc[i,j] = 1 + (float(num_gene.search(l.title).group(1)) -3)
                            break
    return pd_revise
    
    
def find_AMR_detail(accession_name_No_genes, reagent, number, temp):
    '''
    for identification of the antimicrobail resistance genes precisely
    The way of precise identification is to identify genes using the database,
    in which each gene has been modified and numbered by the author of BeMAp. 
    '''
    
    accession = accession_name_No_genes['accession']
    No_genes = accession_name_No_genes['No_genes']
    
    num_gene = re.compile(r'\|\|(.*)?')

    # find AMR genes. cells in no gene is to be 0.01
    detect_AMR = pd.DataFrame(np.array([[0.01]*len(accession) for i in range(max(No_genes))]),columns = accession)
    for i in accession_name_No_genes.index:
        for k in range(accession_name_No_genes.loc[i, 'No_genes']):
            cline = Applications.NcbiblastnCommandline(query = temp + '/each_fasta/' + accession_name_No_genes.loc[i, 'accession'] +'/' + str(k) + '.fasta', db='database/ResFinder_db_detail/' + reagent + '_detail.fsa', outfmt=5)()[0]
            blast_result=NCBIXML.read(io.StringIO(cline))
            description = blast_result.descriptions
            if description != []:
                for l in description:
                    if l.e == 0:
                        detect_AMR.iloc[k,i] = float(num_gene.search(l.title).group(1))
                        break
                    else:
                        detect_AMR.iloc[k,i] = 0
            else:
                detect_AMR.iloc[k,i] = 0

    return detect_AMR
