import logging
import numpy as np
import pandas as pd

fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
logger = logging.getLogger('LogBeMAp').getChild('sub')
logging.basicConfig(level=logging.INFO, format=fmt)

def divide_max_upper(de, gene): # de must be pd.DataFrame and gene must be list
    ranking = pd.DataFrame(np.zeros([len(gene),2]), columns =['max','location'],index=gene, dtype = 'int64')
    for i in gene:
        ranking.loc[i,'max']= (de == i).sum().max()
        ranking.loc[i,'location'] = (de == i).sum().idxmax()
    max_index = ranking['max'].idxmax()
    new_de_1 = de.reindex((de == max_index).sort_values(by = ranking['location'].loc[max_index], ascending=False).index[:ranking['max'].max()])
    del new_de_1[ranking['location'].loc[max_index]]
    return new_de_1

def divide_max_lower(de, gene): # de must be pd.DataFrame and gene must be list
    ranking = pd.DataFrame(np.zeros([len(gene),2]), columns =['max','location'],index=gene, dtype = 'int64')
    for i in gene:
        ranking.loc[i,'max']= (de == i).sum().max()
        ranking.loc[i,'location'] = (de == i).sum().idxmax()
    max_index = ranking['max'].idxmax()
    new_de_2 = de.reindex((de == max_index).sort_values(by = ranking['location'].loc[max_index], ascending=False).index[ranking['max'].max():])
    return new_de_2
    
def check_others_filled(target_centered):
    gene = [x for x in np.unique(target_centered) if np.isnan(x) == False][1:]
    
    sort_by_length = target_centered.isnull().sum(axis='index').sort_values().index
    amr = target_centered.T.reindex(sort_by_length)
    contain_all_index = pd.DataFrame([amr.iloc[0,:]]*len(amr.index),index = amr.index, columns = amr.columns)
    if ((contain_all_index[contain_all_index > 0] == amr).sum() == len(amr.index)).sum() > 1:
        return True
    else:
    	return False
    	
def del_filled_lane(target_centered):
    gene = [x for x in np.unique(target_centered) if np.isnan(x) == False][1:]
    
    sort_by_length = target_centered.isnull().sum(axis='index').sort_values().index
    amr = target_centered.T.reindex(sort_by_length)
    contain_all_index = pd.DataFrame([amr.iloc[0,:]]*len(amr.index),index = amr.index, columns = amr.columns)
    filled_lane = contain_all_index.T[((contain_all_index[contain_all_index > 0] == amr).sum() == len(amr.index))].index
    check_lane = list(target_centered.index)
    for i in filled_lane[:-1]:
        check_lane.remove(i)
    return target_centered.reindex(check_lane)

def align_target_centered(target_centered):
    gene = [x for x in np.unique(target_centered) if np.isnan(x) == False][1:]
    
    sort_by_length = target_centered.isnull().sum(axis='index').sort_values().index
    amr = target_centered.T.reindex(sort_by_length)
    
    du = [[amr]]
    df = [[]]
    k = 1
    while len(du)<len(amr):
        pre_df = du[k-1]
        du.append([])
        du[k] = []
        df.append([])
        df[k] = []
        print(k,'here')
        print(df)
        for i in du[k-1]:
            if len(i) == 0:
                continue
            elif len(i) == 1:
                du[k].append(i)
                df[k].append(list(i.index))
            else:
                cop = i.copy()
                Upper = divide_max_upper(cop,gene)
                du[k].append(Upper)
                df[k].append(list(Upper.index))
                Lower = divide_max_lower(cop,gene)
                if Lower.empty:
                    continue
                else:
                    du[k].append(Lower)
                    df[k].append(list(Lower.index))
        du.append(du[k])
        df.append(df[k])
        after_df = df[k]
        if k != 1:
            if len(pre_df) == len(after_df):
                break
            else:
                k += 1
        else:
            k += 1
    logging.info('Division is stopped at ' + str(k-2))
    
    # collect indexes
    df1 = df
    dfx = []
    
    for i in df1:
        df1_each = []
        for j in i:
            if j != []:
                df1_each.append(j)
        dfx.append(df1_each)
    
    pd_dfx = pd.DataFrame(index= dfx[1][0])
    
    ddfx = dfx.copy()
    del ddfx[0]
    
    for i in range(len(ddfx)):
        for j in range(len(ddfx[i])):
            for l in ddfx[i][j]:
                pd_dfx.loc[l,i] = j
                
    wah = []
    for i in df1[-1]:
        if i != []:
            for j in i:
                wah.append(j)
    new_acc = amr.reindex(wah)
    return new_acc, pd_dfx.sort_values(pd_dfx.columns[-1]).T.reindex(range(k-2)).T