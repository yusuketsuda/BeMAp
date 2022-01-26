from Bio import Phylo
import pandas as pd
import numpy as np
import re
from io import StringIO
from PIL import Image, ImageOps
import matplotlib.pyplot as plt

def make_dendrogram(pd_dfx, tem_store):
    #import a csv file concerning with numbering of each plasmid
    du = pd.DataFrame(pd_dfx.values,index=pd_dfx.index)
    du_sort = du.sort_values(len(pd_dfx.columns)-1)
    
    # unite from back
    # select unique number and, if same, store into dictonary and accession names into key of dictionary. Repeat.
    da = du_sort.loc[:,:len(pd_dfx.columns)]
    for i in reversed(range(1,len(pd_dfx.columns))):
        div_dict_da = {}
        for j in da[i].unique():
            if len(da[da[i] == j].index) == 1:
                div_dict_da[str(list(da[da[i] == j].index)).replace("\\","")] = da.loc[da[i] == j,:i-1].values[0] # remove \. (if not, speed down)
            else:
                div_dict_da[str(list(da[da[i] == j].index)).replace("\\","")] = da.loc[da[i] == j,:i-1].values[0]
        da = pd.DataFrame.from_dict(div_dict_da).T
    all_index  = da.index
    new_index = []
    
    for i in all_index:
        new_index.append(i.translate(str.maketrans({'[':'(',']':')','\'':'','\"':''}))) # convert to newick style
    newick_all = str(new_index).translate(str.maketrans({'[':'(',']':')','\'':''})) # convert to newick style
    
    newick_blank1 = re.sub('[A-Z]*\_', ' ', newick_all)
    newick_blank = re.sub('[A-Z]*[0-9]*\.[0-9]*', ' ', newick_blank1)
    tree = Phylo.read(StringIO(newick_blank), "newick")
    tree.rooted = True
    fig = plt.figure(figsize=(2,100),dpi=100)
    plt.rcParams["font.size"] = 10
    axes = fig.add_subplot(1, 1, 1)
    
    fig.patch.set_facecolor('white') 
    axes.axis('off')
    
    Phylo.draw(tree,axes=axes,do_show = False)
    plt.savefig( tem_store + 'pre_de.png',bbox_inches='tight', pad_inches=0)
    
    im = Image.open( tem_store + 'pre_de.png')
    change = ImageOps.flip(im.rotate(90,expand=True))
    
    pd_im = pd.DataFrame(np.array(change.convert('L')))
    de = pd_im[pd_im == 0].dropna(how='all',axis=0).dropna(how='all',axis=1)
    mar = 20
    da = pd_im.loc[de.index[0]-3:de.index[-1]+1,de.columns[0]-mar:de.columns[-1]+mar]
    
    save_img = Image.fromarray(da.fillna(255).to_numpy().astype('uint8'))
    save_img.save('figure/dendrogram.png')
    
    plt.clf()
    plt.close()
