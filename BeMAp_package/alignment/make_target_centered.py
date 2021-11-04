from Bio import SeqIO
import pandas as pd
import numpy as np
import math
import re
import os 

def invert_new(pdDataFrame, No_genes, pd_column_index, gene_description):
    pdDataFrame_dict = pdDataFrame.to_dict()
    new_pdDataFrame_dict = {}
    for i in No_genes.index:
        new_pdDataFrame_dict[i] = list(pdDataFrame_dict[i].values())[:No_genes[i]]
    
    complement = re.compile(r'strand_(.*)?\|\|')
    
    pdDataFrame_dict_inverted = {}
    for i in pd_column_index.columns:
        for j in pd_column_index.index:
            if 1 <= pd_column_index.at[j,i] < 2:
                    if complement.search(gene_description.at[j,i]).group(1) == '-1': 
                        pdDataFrame_dict_inverted[i] = list(new_pdDataFrame_dict[i])[::-1]
                    else:
                        pdDataFrame_dict_inverted[i] = list(new_pdDataFrame_dict[i])
            else:
                continue
    return pd.DataFrame.from_dict(pdDataFrame_dict_inverted, orient='index').T

def make_target_centered(list_AMR, temp):
    AMR = list_AMR[0].copy()
    for i in range(len(list_AMR))[::-1]:
        AMR.update(list_AMR[i][list_AMR[i].astype(int) == list_AMR[i].astype(int).max().max()])
    AMR.replace({0.01:np.nan}, inplace=True)

    No_genes = (~(AMR.isnull())).sum()
    
    gene_dict = {}
    description = {}
    gene_name = re.compile(r'\|\|\s(.*)?')
    for i in No_genes.index:
        each_gene = []
        each_description = []
        for j in range(No_genes[i]):
            record = SeqIO.read(temp + '/each_fasta/' + i +'/' + str(j) + '.fasta','fasta')
            each_description.append(record.description)
            each_gene.append(gene_name.search(record.description).group(1))
        gene_dict[i] = each_gene
        description[i] = each_description
    gene = pd.DataFrame.from_dict(gene_dict, orient='index').T
    gene_description = pd.DataFrame.from_dict(description, orient='index').T
    
    AMR_inverted = invert_new(AMR, No_genes, AMR, gene_description)
    
    target_len = pd.DataFrame(columns=['complete sequence','position','len','seq'])
    AMR_inverted_int = AMR_inverted.fillna(0).astype(int)
    for i in AMR_inverted_int.T.index:
        position = []
        length = []
        sequence = []
        for j in AMR_inverted_int[i][AMR_inverted_int[i]==1].index:
            record = SeqIO.read(temp + '/each_fasta/' + i +'/' + str(j) + '.fasta','fasta')
            position.append(j)
            length.append(len(record.seq))
            sequence.append(str(record.seq))
        target_len.loc[i,'position'] = position
        target_len.loc[i,'len'] = length
        target_len.loc[i,'seq'] = sequence
    
    for i in No_genes.index:
        record = SeqIO.read(temp + '/fasta/'+ i + '.fasta','fasta')
        if re.compile(r'complete sequence').search(record.description) or re.compile(r'complete genome').search(record.description): 
            if not re.compile(r'cds').search(record.description):
                target_len.loc[i,'complete sequence'] = 'yes'
    
    # targetを先頭にする
    first = {}
    for i in AMR_inverted.columns:
        position = target_len.loc[i,'position'][0]
        upper = list(AMR_inverted.T.loc[i,position:].dropna())
        lower = list(AMR_inverted.T.loc[i,:position-1].dropna())
        first[i] = upper + lower
    target_first = pd.DataFrame.from_dict(first,orient='index')
    
    gene_inverted = invert_new(gene, No_genes, AMR, gene_description)
    
    # geneの種類についてtargetを先頭にする
    first_gene = {}
    for i in AMR_inverted.columns:
        position = target_len.loc[i,'position'][0]
        upper = list(gene_inverted.T.loc[i,position:].dropna())
        lower = list(gene_inverted.T.loc[i,:position-1].dropna())
        first_gene[i] = upper + lower
    gene_target_first = pd.DataFrame.from_dict(first_gene,orient='index')
    
    # 最大のレーン決める
    mid = pd.DataFrame(columns=['mid'])
    for i in target_first.index:
        mid.loc[i,'mid'] = math.ceil((len(target_first.loc[i,:]) - target_first.loc[i,:].isnull().sum())/2)
    maximum = mid.max()
    
    # とりあえずplasmidを半分にして、targetをだいたい中央に持ってくる
    center = {}
    for i in target_first.index:
        Mid = int(mid.loc[i,'mid'])
        center[i] = [None] * (int(maximum[0])-Mid) + list(target_first.loc[i,Mid+1:].dropna()) + list(target_first.loc[i,:Mid])
    modified_AMR = pd.DataFrame.from_dict(center, orient='index').T
    
    # とりあえずplasmidを半分にして、targetをだいたい中央に持ってくる
    gene_center = {}
    for i in gene_target_first.index:
        Mid = int(mid.loc[i,'mid'])
        gene_center[i] = [None] * (int(maximum[0])-Mid) + list(gene_target_first.loc[i,Mid+1:].dropna()) + list(gene_target_first.loc[i,:Mid])
    modified_gene = pd.DataFrame.from_dict(gene_center, orient='index').T
    
    # targetが中央になるようにする。これで完成
    # completeでないdata setはそのままの順序で挿入しておく。
    if len(modified_AMR.index) > 10:
    	around_max_target10 = modified_AMR.loc[int(maximum[0])-5:int(maximum[0])+5,:].fillna(0).astype(int)
    else:
    	around_max_target10 = modified_AMR.fillna(0).astype(int)
    count_target = (around_max_target10==1).sum(axis='columns')
    count_target_max = count_target.max()
    for i in list(count_target[count_target == count_target_max].index):
    	target_max_row = i
    if count_target[target_max_row+1] > 0:   # 1つ下の行に1があれば、そちらに合わせる。
        target_max_row += 1

    new_modified = {}
    gene_new_modified = {}
    for i in modified_AMR.columns:
        if target_len.loc[i,'complete sequence'] == 'yes':
            align_AMR = list(modified_AMR.T.loc[i,:])
            align_gene = list(modified_gene.T.loc[i,:])
            d = 0
            while modified_AMR.fillna(0).astype(int).T.loc[i,target_max_row-d] != 1:
            	d += 1
            new_modified[i] = [np.nan] * d + align_AMR
            gene_new_modified[i] = [np.nan] * d + align_gene
        else:
            align_AMR = list(AMR_inverted.T.loc[i,:])
            align_gene = list(gene_inverted.T.loc[i,:])
            new_modified[i] = [np.nan] * (target_max_row - target_len.loc[i,'position'][0]) + align_AMR
            gene_new_modified[i] = [np.nan] * (target_max_row - target_len.loc[i,'position'][0]) + align_gene
    pd_new_modified = pd.DataFrame.from_dict(new_modified, orient='index').T
    pd_gene_new_modified = pd.DataFrame.from_dict(gene_new_modified, orient='index').T
    
    return pd_new_modified, pd_gene_new_modified
