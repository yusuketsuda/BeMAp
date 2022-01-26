#!/usr/bin/python

import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger('LogBeMAp').getChild('sub')
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
logging.basicConfig(level=logging.INFO, format=fmt)

def make_each_fasta(accession, genbank, fasta, temp,s=True): 
    '''
    s : If True, save fasta files (default: True)
    
    In our thesis, all datasets were annotated by prokka, so the style of genbank were consistent, leading good analysis and mapping.    
    In fact, you do not necessarily annotate all datasets and you can use datasets downloaded from websites or made by yourself.
    However, some of such datasets may contain non-usual features so that some of genes cannot be analyzed.
    
    '''
    try:
        record = SeqIO.read(genbank, 'genbank')
        record_fasta = SeqIO.read(fasta,'fasta')
        product_gene_note = []
        location = []
        strand = []
        for j in range(len(record.features)):
            location_each = []
            if record.features[j].type == 'gene':
                if j != range(len(record.features))[-1]:
                    if record.features[j].location != record.features[j+1].location:
                        if record.features[j].type == 'gene':
                            if 'product' in record.features[j].qualifiers:
                                product_gene_note.append(record.features[j].qualifiers['product'][0])    
                                strand.append(record.features[j].location.strand)
                            elif 'gene' in record.features[j].qualifiers:
                                product_gene_note.append(record.features[j].qualifiers['gene'][0])
                                strand.append(record.features[j].location.strand)
                            elif 'note' in record.features[j].qualifiers:
                                product_gene_note.append(record.features[j].qualifiers['note'][0])
                                strand.append(record.features[j].location.strand)
                            elif 'locus_tag' in record.features[j].qualifiers:
                                product_gene_note.append(record.features[j].qualifiers['locus_tag'][0])
                                strand.append(record.features[j].location.strand)
                            else:
                                print('gene',record.features[j].qualifiers)
                        for l in range(len(record.features[j].location.parts)):
                            if 'product' in record.features[j].qualifiers:
                                location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                            elif 'gene' in record.features[j].qualifiers:
                                location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                            elif 'note' in record.features[j].qualifiers:
                                location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                            elif 'locus_tag' in record.features[j].qualifiers:
                                location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                            else:
                                print('gene_location',record.features[j].qualifiers)
                        location.append(location_each)
            elif record.features[j].type == 'CDS':
                    if 'product' in record.features[j].qualifiers:
                        product_gene_note.append(record.features[j].qualifiers['product'][0])    
                        strand.append(record.features[j].location.strand)
                    elif 'gene' in record.features[j].qualifiers:
                        product_gene_note.append(record.features[j].qualifiers['gene'][0])
                        strand.append(record.features[j].location.strand)
                    elif 'note' in record.features[j].qualifiers:
                        product_gene_note.append(record.features[j].qualifiers['note'][0])
                        strand.append(record.features[j].location.strand)
                    else:
                        print(record.features[j].qualifiers)
                    for l in range(len(record.features[j].location.parts)):
                        if 'product' in record.features[j].qualifiers:
                            location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                        elif 'gene' in record.features[j].qualifiers:
                            location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                        elif 'note' in record.features[j].qualifiers:
                            location_each.append([record.features[j].location.parts[l].start,record.features[j].location.parts[l].end])
                        else:
                            print('location',record.features[j].qualifiers) 
                    location.append(location_each)
            else:
                continue

        try:
            product_location = pd.DataFrame({'product':product_gene_note, 'location':location,'strand':strand})
            for k in product_location.index:
                if len(product_location.loc[k,'location']) == 1:
                    if s:
                        seq =SeqRecord(record_fasta[product_location.loc[k,'location'][0][0]:product_location.loc[k,'location'][0][1]].seq, id = accession + '_CDS_'+ str(k) + '_strand_'+str(product_location.loc[k,'strand'])+ '||', description = product_location.loc[k,'product'])
                        SeqIO.write(seq, temp + '/each_fasta/' + accession + '/' + str(k) +'.fasta','fasta')
                else:
                    if s:
                        seq = SeqRecord(record_fasta[product_location.loc[k,'location'][0][0]:product_location.loc[k,'location'][0][1]].seq+record_fasta[product_location.loc[k,'location'][1][0]:product_location.loc[k,'location'][1][1]].seq, id = accession + '_CDS_'+ str(k)+'_strand_'+str(product_location.loc[k,'strand'])+ '||', description = product_location.loc[k,'product'])
                        SeqIO.write(seq, temp + '/each_fasta/' + accession + '/' + str(k) +'.fasta','fasta')
            return [accession, k+1]
        except:
            logging.warning(accession + ' is excluded')
    except:
        logging.warning(accession + ' has no record')
