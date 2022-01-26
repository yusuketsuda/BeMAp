import pandas as pd
import numpy as np
import io
import os
import re
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import Applications
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger('LogBeMAp').getChild('sub')
fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
logging.basicConfig(level=logging.INFO, format=fmt)

def del_kakko(name):
    if '(' in name:
        to_Series = pd.Series(list(name))
        position_kakko = to_Series[to_Series == '('].index
        return ''.join(list(to_Series[:position_kakko[0]]))
    else:
        return name

def Inc_extract(name):
    if name == 'IncB/O/K/Z':
        return 'IncB/O/K/Z'
    elif name == 'FII':
        return 'IncF'
    elif name == 'FIA':
        return 'IncF'
    else:
        return ''.join(list(name)[:4])

def only_inc(name):
    under = re.compile(r'(.*?)\_')
    if 'p' == list(name)[0] or 'C' == list(name)[0] or 'r' == list(name)[0]:
        return 
    else:
        return under.search(name).group(1)
    
def summarize_inc(list_inc):
    IncA_C = {'IncA','IncC'}
    IncL_M = {'IncL','IncM'}
    IncHI  = {'IncH'}
    
    if len(list_inc) > 0:
        each_inc = []
        for i in list_inc:
            if only_inc(i):
                each_inc.append(Inc_extract(del_kakko(i)))
        set_inc_group = set(each_inc)
        
        if 'IncA' in set_inc_group and 'IncC' in set_inc_group:
            set_inc_group = set_inc_group - {'IncA','IncC'} | {'IncA/C'}
        if 'IncL' in set_inc_group and 'IncM' in set_inc_group:
            set_inc_group = set_inc_group - {'IncL','IncM'} |{'IncL/M'}
        if 'IncH' in set_inc_group:
            set_inc_group = set_inc_group - {'IncH'} |{'IncHI'}
            
        if len(list(set_inc_group)) == 1:
            return list(set_inc_group)[0]
        else:
            return ','.join(list(set_inc_group))

# select dataset typed as 'plasmid'
def find_plasmid(pre_accession, genbank_files, temp):
    all_plasmid_accession = []
    all_plasmid_genbank = []
    for i in range(len(pre_accession)):
        record = SeqIO.read(temp + '/genbank/'+ genbank_files[i],'genbank')
        features = record.features
        if 'plasmid' in features[0].qualifiers:
            all_plasmid_accession.append(pre_accession[i])
            all_plasmid_genbank.append(genbank_files[i])
        else:
            logging.warning(pre_accession[i] + ' was excluded: non-\'plasmid\'')
    
    # exclude datasets without information of fragments
    plasmid_accession = []
    plasmid_genbank = []
    for i in range(len(all_plasmid_accession)):
        record = SeqIO.read(temp + '/genbank/'+ all_plasmid_genbank[i],'genbank')
        features = record.features
        if len(features) > 1:
            plasmid_accession.append(all_plasmid_accession[i])
            plasmid_genbank.append(all_plasmid_genbank[i])
        else:
            logging.warning(all_plasmid_accession[i] + ' was excluded: no information')
    return plasmid_accession, plasmid_genbank
            
def find_data_containing_target(pd_target):
    accession_containing_target = pd_target.columns[pd_target[pd_target == 1].sum() > 0]
    return list(accession_containing_target)

def summarize_plasmids(accession, loc_genbank, result_plasmid, result_contain_target, temp):
    inc = re.compile(r'\|\d*\s(.*)?')    
    summary = pd.DataFrame(index=accession, columns=['type','name','size(bp)','organism','host','country','position of target','include?','Inc group','Inc group (Results by PlasmidFinder)'])
    plasmid_include_target = list(set(result_plasmid) & set(result_contain_target))
    x = 0
    for i in accession:
        genbank = temp + '/genbank/'+ loc_genbank[x]
        fasta   = temp + '/fasta/' + i + '.fasta'
        record_gb = SeqIO.read(genbank,'genbank')
        record_fa = SeqIO.read(fasta,'fasta')
        # include?
        if record_gb.id in plasmid_include_target:
            summary.loc[i, 'include?'] = 'yes'
        else:
            summary.loc[i, 'include?'] = 'no'    
        # size
        summary.loc[i, 'size(bp)'] = len(record_fa)
        # type
        if 'plasmid' in record_gb.features[0].qualifiers:
            summary.loc[i, 'type'] = 'plasmid'
            summary.loc[i, 'name'] = record_gb.features[0].qualifiers['plasmid'][0]
        # country
        if 'country' in record_gb.features[0].qualifiers:
            summary.loc[i, 'country'] = record_gb.features[0].qualifiers['country'][0]
        # host
        if 'host' in record_gb.features[0].qualifiers:
            summary.loc[i, 'host'] = record_gb.features[0].qualifiers['host'][0]    
        # organism
        if 'organism' in record_gb.features[0].qualifiers:
            summary.loc[i,'organism'] = record_gb.features[0].qualifiers['organism'][0]
        # position of target
        cline = Applications.NcbiblastnCommandline(query=fasta, db='database/target/target.fsa',evalue = 0.1,outfmt=5)()[0]
        blast_result =NCBIXML.read(io.StringIO(cline))
        descriptions = blast_result.descriptions
        alignments = blast_result.alignments
        position = []
        for alignment in alignments:
            for hsp in alignment.hsps:
                position.append(str(hsp.query_start) + ':' + str(hsp.query_end))
        if len(position) > 1:
            unite_position = position[0]
            for p in range(len(position)-1):
                unite_position += ', ' + position[p+1] 
            summary.loc[i, 'position of target'] = unite_position
        elif len(position) == 1:
            summary.loc[i, 'position of target'] = position[0]
        # classfication of Inc groups
        if 'plasmid' in record_gb.features[0].qualifiers:
            inc_each = []
            cline = Applications.NcbiblastnCommandline(query=fasta, db='database/PlasmidFinder_db/enterobacteriaceae.fsa',evalue = 0.001,perc_identity=90,outfmt=5)()[0]
            blast_result =NCBIXML.read(io.StringIO(cline))
            descriptions = blast_result.descriptions
            for des in descriptions:
                inc_each.append(inc.search(des.title).group(1))
            sorted_inc_each = sorted(inc_each)
            summary.loc[i, 'Inc group'] = summarize_inc(sorted_inc_each)
            if len(sorted_inc_each) > 1:
                new_inc_each = sorted_inc_each[0]
                for k in range(len(sorted_inc_each)-1):
                    new_inc_each += ', ' + str(sorted_inc_each[k+1]) 
                summary.loc[i, 'Inc group (Results by PlasmidFinder)'] = new_inc_each
            elif len(sorted_inc_each) == 1:
                summary.loc[i, 'Inc group (Results by PlasmidFinder)'] = sorted_inc_each[0]
        x += 1
    return summary
