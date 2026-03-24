import os
import sys
import pandas as pd  
from scipy.sparse import csr_matrix  
import numpy as np  
import matplotlib.pyplot as plt
import seaborn as sns
#import dask.dataframe as dd  
sys.path.append("~/NPC_project/BMK/pipeline/amplicon")
import preprocess_utils as pt
wd="~/NPC_project/BMK/tissue/brain/P5/S3000/amplicon_all/tmp" 
wl=pd.read_table("~/NPC_project/P045_whitelist/P045/whitelist.tsv",header=None)
input=os.path.join(wd,"collapsed_clonebc_umi_cell.txt")
cr_df=pd.read_table(input,header=None)[[4,1,0,2,3]]
cr_df.columns = ['umi','clonebc','cell','cell_x','cell_y']
cr_df[:5]

#filter with whitelist
df0=cr_df[cr_df['clonebc'].isin(wl[0])]
cr_df.shape,df0.shape

print(str(df0.shape[0])+"/"+str(cr_df.shape[0])+" cell-clonebc are retained after filtering with whitelist")
df1=pt.count_umi(df0,['cell','clonebc','cell_x','cell_y'
                     ],stat='uniq')
pt.get_scatter_plot(df1,group_col='umi')

df1a = df1#[df1['clonebc'].isin(dominant_clonebc)]
print("celll-clonebc has "+str(df1a.shape[0]))
pt.get_scatter_plot(df1a,group_col='umi')
df1a.to_csv(os.path.join(wd,"qc1_cell_clonebc_umi.tsv"),sep='\t')
print("file is saved to "+wd+"/qc1_cell_clonebc_umi.tsv")

df1b = df1a#[(df1a['umi'] > 1)]
print("celll-clonebc has "+str(df1b.shape[0]))
pt.get_scatter_plot(df1b,group_col='umi')
#df1b.to_csv(os.path.join(wd,"qc2_cell_clonebc_umi.tsv"),sep='\t')
print("file is saved to "+wd+"/qc2_cell_clonebc_umi.tsv")

df1c=pt.clonebc_umi_percent_percell(df1b,thr_low=0.1,thr_high=1.1,umi_thr_low=0,umi_thr_high=2**8)
pt.get_scatter_plot(df1c[df1c['confident_clonebc']],group_col='umi')

df2=pt.cell_conf_clonebc_percent(df1c,thr=0.5)
df3=pt.clonebc_conf_cell_percent(df2,thr=0.5)
final_df=df3[df3['confident_clonebc']  & df3['confident_cell']  ]#& df3['confident_clonebc_by_truecell_percent']

print("final confident celll-clonebc has "+str(final_df.shape[0]))
pt.get_scatter_plot(final_df,group_col='umi')

final_df.to_csv(os.path.join(wd,"qc3_cell_clonebc_umi.tsv"),sep='\t')
print("file is saved to "+wd+"/qc3_cell_clonebc_umi.tsv")

### read, cell, clonebc statistics
clonebc_stat=pd.read_table(os.path.join(wd,"clonebc.stat.txt"),index_col=0#,header=None
                          )
cell_stat=pd.read_table(os.path.join(wd,"cell.stat.txt"),index_col=0#,header=None
                       )
print(clonebc_stat)
print(cell_stat)

ncell1=len(np.unique(df0[['cell']]))
nbc1=len(np.unique(df0[['clonebc']]))
ncell2=len(np.unique(final_df[['cell']]))
nbc2=len(np.unique(final_df[['clonebc']]))
ncell1,nbc1,ncell2,nbc2

ncell2=len(np.unique(final_df[['cell']]))
nbc2=len(np.unique(final_df[['clonebc']]))
ncell1,nbc1,ncell2,nbc2

clonebc_stat.loc[4] = ['whitelist', nbc1]
clonebc_stat.loc[5] = ['confident', nbc2]
clonebc_stat

cell_stat.loc[2] = ['whitelist', ncell1]
cell_stat.loc[3] = ['confident', ncell2]
cell_stat

clonebc_stat.to_csv(os.path.join(wd,'clonebc.stat.txt'),sep='\t')
cell_stat.to_csv(os.path.join(wd,'cell.stat.txt'),sep='\t')
