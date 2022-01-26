import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import numpy as np

def make_AMR_legend(color, temp_dir):
    for i in color.index:
        plt.scatter([],[],marker='s',s=300,c= '#' + color.loc[i,1])
    plt.axis('off')
    plt.rcParams["legend.loc"] = "lower left" 
    plt.rcParams["font.size"] = 20
    plt.rcParams["legend.fancybox"] = False
    plt.rcParams["legend.labelspacing"] = 0
    plt.rcParams["legend.handletextpad"] = -0.2
    legend = plt.legend(list(color[0]),title='Types of the antimicrobial resistance genes',ncol=14)
    legend._legend_box.align = 'left'
    plt.savefig(temp_dir + 'pre_AMR.png',bbox_inches='tight',facecolor='white')
    plt.clf()
    plt.close()
    
    np_im = np.array(Image.open(temp_dir + 'pre_AMR.png'))
    
    pd_im1 = pd.DataFrame(np_im[:,:,0])
    pd_im2 = pd.DataFrame(np_im[:,:,1])
    pd_im3 = pd.DataFrame(np_im[:,:,2])
    pd_im4 = pd.DataFrame(np_im[:,:,3])
    de = pd_im1[pd_im1 != 255].dropna(how='all',axis=0).dropna(how='all',axis=1)
    
    pd_to_np = np.array(
                [pd_im1.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im2.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im3.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im4.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T]
                )
    
    save_img = Image.fromarray(pd_to_np.astype('uint8').T)
    save_img.save('figure/legend/AMR.png')
    
    plt.close()
    

def make_property_legend(prop, prop_color, temp):
    leg = []
    for i in range(len(prop_color)):
        for j in range(4):
            if j%4 == 0:
                plt.scatter([],[],marker='s',alpha=0)
                leg.append('')
            else:
                if j%4 == 2:
                    leg.append(prop_color.iloc[i,0])
                    plt.scatter([],[],marker='s',s=100,c='#' + prop_color.iloc[i,1])
                else:
                    leg.append('')
                    plt.scatter([],[],marker='s',s=100,c='#' + prop_color.iloc[i,1])
    plt.axis('off')
    plt.rcParams["legend.labelspacing"] = -0.3
    plt.rcParams["font.size"] = 20
    plt.rcParams["legend.handletextpad"] = -0.2
    legend = plt.legend(leg,title='Types of ' + prop, title_fontsize='medium',ncol=len(prop_color))
    legend._legend_box.align = 'left'
    plt.savefig(temp + '/'+ prop + '.png',bbox_inches='tight',facecolor='white')
    plt.clf()
    plt.close()
    
    np_im = np.array(Image.open(temp + '/' + prop + '.png'))
    
    pd_im1 = pd.DataFrame(np_im[:,:,0])
    pd_im2 = pd.DataFrame(np_im[:,:,1])
    pd_im3 = pd.DataFrame(np_im[:,:,2])
    pd_im4 = pd.DataFrame(np_im[:,:,3])
    de = pd_im1[pd_im1 != 255].dropna(how='all',axis=0).dropna(how='all',axis=1)
    
    pd_to_np = np.array(
                [pd_im1.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im2.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im3.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T,
                 pd_im4.loc[de.index[0]:de.index[-1],de.columns[0]:de.columns[-1]].fillna(255).to_numpy().T]
                )
    
    save_img = Image.fromarray(pd_to_np.astype('uint8').T)
    save_img.save('figure/legend/' + prop + '.png')
    
    plt.close()
