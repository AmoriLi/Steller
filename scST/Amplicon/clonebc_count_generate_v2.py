import os
import sys
import pandas as pd  
from scipy.sparse import csr_matrix  
import numpy as np  
#input=sys.argv[1] #e.g. "~/NPC_project/BMK/tissue/brain/E15.5/2024.5.7/S3000-A4/amplicon/result/CloneBC_v4/1_preprocess/spbc2.umi2.clonebc"
wd=sys.argv[1] #"~/NPC_project/BMK/tissue/brain/E15.5/2024.5.7/S3000-A4/amplicon/tmp"
input=os.path.join(wd,"umi2.clonebc.cellmeta.nreadcut.tsv")#"umi2_clonebc_cell.txt")#"spbc2.umi2.clonebc.txt")
# 
def count_umi(data,group):  
    umi_counts = data.groupby(group)['umi'].nunique().reset_index()  
    group.append('umi_count')
    umi_counts.columns = group
    return umi_counts 

data = pd.read_table(input,header=None)[[4,1,0,2,3]]#[[5,1,2,3,4]] 
data.columns=['umi','clonebc','cell','cell_x','cell_y'] #need correct

umi_counts = count_umi(data,group=["clonebc"]) 
umi_counts.to_csv(os.path.join(wd,"clonebc_count.txt"),sep='\t',header=None,index=None)