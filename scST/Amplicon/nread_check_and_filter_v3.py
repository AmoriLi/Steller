##2025-10-19
#cel total reads check and filter. >1000
import os
import sys
import pandas as pd  
from scipy.sparse import csr_matrix  
import numpy as np  
import matplotlib.pyplot as plt
#import dask.dataframe as dd  
sys.path.append("~/NPC_project/BMK/pipeline/amplicon")
import preprocess_utils as pt

wd=sys.argv[1] #"~/NPC_project/BMK/tissue/brain/E15.5/2024.5.7/S3000-A4/amplicon/tmp" #sys.argv[1] #
thr1=float(sys.argv[2])
thr2=float(sys.argv[3])
input=os.path.join(wd,"spbc2.umi2.clonebc.cellmeta.tsv")

data = pd.read_table(input,header=None)[[4,15,16,17,18]]
data.columns = ['umi','clonebc','cell','cell_x','cell_y'#,'spbc1','spbc2','spbc3','spbc1_seq','spbc2_seq','spbc3_seq','spbc_seq','seq','umi_seq_group','clonebc','cell','cell_x','cell_y','reads','umi_passed_nread'#,'swapped_ncells','swapped','umi_clonebc_group','swapped_ncells2'
               ]
#data=data[['id','spbc','spbc_umi0','umi0','umi','spbc_seq','seq','umi_seq_group','clonebc','cell','cell_x','cell_y']]
#data[:5]
#transform('size')
df=pt.cell_clonebc_read_check(data,inplace=True,read_lowthr=thr1,read_highthr=2 ** 20,addlog=True,
                                  saveplot=True,outdir=wd)
df2=pt.cell_read_check(df,inplace=True,read_lowthr=thr2,read_highthr=2 ** 20,addlog=True,saveplot=True,outdir=wd)

df2.to_csv(os.path.join(wd,"umi2.clonebc.cellmeta.nreadcut.tsv"),header=None,sep='\t',index=None)

#print(str(df[df['umiread_passed']].shape[0])+" reads are retained after umi_nread cutoff!")