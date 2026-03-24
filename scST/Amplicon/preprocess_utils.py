import os
import sys
import pandas as pd  
from scipy.sparse import csr_matrix  
import numpy as np  
import matplotlib.pyplot as plt
#import dask.dataframe as dd  

def get_hist_threshold_percent_plot(df, col, threshold_low,threshold_high=None,bins=100,addlog=False,save=False,filename='read',outdir=None):
    NREAD1 = threshold_low
    NREAD2 = threshold_high

    plt.figure(figsize=(8, 6))  
    data=df[col]
    if addlog:
        data=np.log2(df[col])  
        NREAD1=np.log2(NREAD1)
        NREAD2=np.log2(NREAD2)

    plt.hist(data, bins=bins, color='grey', edgecolor='black')  
    if NREAD2 is None:
        above_threshold_count = (data > NREAD1).sum()#shape[0]
        total_count = df.shape[0]
        above_threshold_percentage = (above_threshold_count / total_count) * 100
    else:
        above_threshold_count = ((data > NREAD1) & (data < NREAD2)).sum()#shape[0]
        total_count = df.shape[0]
        above_threshold_percentage = (above_threshold_count / total_count) * 100

    if addlog:
        plt.axvline(x=NREAD1, color='red', linestyle='--', label=f'Threshold lower: {2 ** NREAD1}')
        if not NREAD2 is None:
            plt.axvline(x=NREAD2, color='red', linestyle='--', label=f'Threshold upper: {2 ** NREAD2}')
        plt.legend()
        if NREAD2 is None:
            plt.text(0.3, 0.5, f'Above {2 ** NREAD1}: {above_threshold_count} ({above_threshold_percentage:.2f}%)', transform=plt.gca().transAxes)
        else:
            plt.text(0.3, 0.5, f'Above {2 ** NREAD1} and below {2 ** NREAD2}: {above_threshold_count} ({above_threshold_percentage:.2f}%)', transform=plt.gca().transAxes)
    else:
        plt.axvline(x=NREAD1, color='red', linestyle='--', label=f'Threshold lower: {NREAD1}')
        if not NREAD2 is None:
            plt.axvline(x=NREAD2, color='red', linestyle='--', label=f'Threshold upper: {NREAD2}')
        plt.legend()
        if NREAD2 is None:
            plt.text(0.3, 0.5, f'Above {NREAD1}: {above_threshold_count} ({above_threshold_percentage:.2f}%)', transform=plt.gca().transAxes)
        else:
            plt.text(0.3, 0.5, f'Above {NREAD1} and below {NREAD2}: {above_threshold_count} ({above_threshold_percentage:.2f}%)', transform=plt.gca().transAxes)

    plt.title(col+' Distribution\nlog2='+str(addlog))
    plt.xlabel(col+' Values')
    plt.ylabel('Frequency')
    #plt.show()
    if save:
        plt.savefig(os.path.join(outdir,filename+'.hist.stat.png'))
        plt.clf()
    else:
        plt.show()

def umi_read_check(df,inplace,#read_group,
                     read_lowthr,read_highthr,addlog,saveplot,outdir):
    raw_df=df
    bc_umi_nread = df#[['umi']]#.drop_duplicates()
    #bc_umi_nread.columns = ['umi','clonebc','cell']
    #bc_umi_nread['read'] = 1
    #bc_umi_nread[:5]
    bc_umi_nread['reads'] = bc_umi_nread.groupby('umi')['umi'].transform('size')#.agg(reads=(read_group, lambda x: x.sum())).reset_index() #.transform('size')
    bc_umi_nread_uniq = bc_umi_nread.drop_duplicates()
    #bc_umi_nread.shape,bc_umi_nread_uniq.shape
    get_hist_threshold_percent_plot(bc_umi_nread_uniq,'reads',read_lowthr,read_highthr,addlog=addlog,bins=100,save=saveplot,outdir=outdir)
    
    if inplace:
        raw_df = bc_umi_nread
        raw_df['umiread_passed'] = ((raw_df['reads'] > read_lowthr) & (raw_df['reads'] < read_highthr))
        #raw_df[raw_df['umi'].isin(df_sub['umi'])]['qc1_filtered']
        raw_df = raw_df[raw_df['umiread_passed']]
        return raw_df
    
def cell_clonebc_umi_read_check(df,inplace,#read_group,
                     read_lowthr,read_highthr,addlog,saveplot,outdir):
    raw_df = df
    bc_umi_nread = df[['umi','cell','clonebc','cell_x','cell_y']]#.drop_duplicates()
    bc_umi_nread['reads'] = bc_umi_nread.groupby(['umi','clonebc','cell','cell_x','cell_y'])['umi'].transform('size')#.agg(reads=(read_group, lambda x: x.sum())).reset_index() #.transform('size')
    bc_umi_nread_uniq = bc_umi_nread.drop_duplicates()
    #bc_umi_nread.shape,bc_umi_nread_uniq.shape
    get_hist_threshold_percent_plot(bc_umi_nread_uniq,'reads',read_lowthr,read_highthr,addlog=addlog,bins=100,save=saveplot,outdir=outdir)
    if inplace:
        raw_df = bc_umi_nread
        raw_df['read_passed'] = ((raw_df['reads'] > read_lowthr) & (raw_df['reads'] < read_highthr))
        raw_df = raw_df[raw_df['read_passed']]#[raw_df['umi'].isin(df_sub['umi'])]['qc1_filtered']
        return raw_df

def cell_clonebc_read_check(df,inplace,#read_group,
                     read_lowthr,read_highthr,addlog,saveplot,outdir):
    raw_df = df
    bc_umi_nread = df[['cell','clonebc','cell_x','cell_y','umi']]#.drop_duplicates()
    
    bc_umi_nread['reads'] = bc_umi_nread.groupby(['clonebc','cell','cell_x','cell_y'])[['umi']].transform('size')#.agg(reads=(read_group, lambda x: x.sum())).reset_index() #.transform('size')
     
    bc_umi_nread_uniq = bc_umi_nread.drop_duplicates()
    #bc_umi_nread.shape,bc_umi_nread_uniq.shape
    get_hist_threshold_percent_plot(bc_umi_nread_uniq,'reads',read_lowthr,read_highthr,addlog=addlog,bins=100,save=saveplot,filename='cell_clonebc_read',outdir=outdir)
    
    if inplace:
        #nread_to_umi = pd.Series(bc_umi_nread_uniq['reads'].values, index=bc_umi_nread_uniq['umi']).to_dict()
        #raw_df['umi_nread']=raw_df['umi'].map(nread_to_umi)
        raw_df = bc_umi_nread
        raw_df['read_passed'] = ((raw_df['reads'] > read_lowthr) & (raw_df['reads'] < read_highthr))
        raw_df = raw_df[raw_df['read_passed']]#[raw_df['umi'].isin(df_sub['umi'])]['qc1_filtered']
        return raw_df 
    
def cell_read_check(df,inplace,#read_group,
                     read_lowthr,read_highthr,addlog,saveplot,outdir):
    raw_df = df
    bc_umi_nread = df[['cell','clonebc','cell_x','cell_y','umi']]#.drop_duplicates()
    
    bc_umi_nread['reads'] = bc_umi_nread.groupby(['cell','cell_x','cell_y'])[['clonebc','umi']].transform('size')#.agg(reads=(read_group, lambda x: x.sum())).reset_index() #.transform('size')
    bc_umi_nread_uniq = bc_umi_nread.drop_duplicates()
    #bc_umi_nread.shape,bc_umi_nread_uniq.shape
    get_hist_threshold_percent_plot(bc_umi_nread_uniq,'reads',read_lowthr,read_highthr,addlog=addlog,bins=100,save=saveplot,filename='cell_read',outdir=outdir)

    if inplace:
        #nread_to_umi = pd.Series(bc_umi_nread_uniq['reads'].values, index=bc_umi_nread_uniq['umi']).to_dict()
        #raw_df['umi_nread']=raw_df['umi'].map(nread_to_umi)
        raw_df = bc_umi_nread
        raw_df['read_passed'] = ((raw_df['reads'] > read_lowthr) & (raw_df['reads'] < read_highthr))
        raw_df = raw_df[raw_df['read_passed']]#[raw_df['umi'].isin(df_sub['umi'])]['qc1_filtered']
        return raw_df
    
def get_scatter_plot(df,x='cell_x',y='cell_y',group_col='conf_clonebc_percent',saveplot=False,outdir=None):
    plt.figure(figsize=(8, 6)) 
    plt.scatter(df[x], df[y],s=0.5, c=df[group_col], cmap='viridis')

    plt.colorbar()

    plt.title(group_col)
    plt.xlabel('x_axis')
    plt.ylabel('y_axis')

    if saveplot:
        plt.savefig(os.path.join(outdir,'cell_conf_clonebc_percent.spatial.png'))
        plt.clf()
    else:
        plt.show()

def count_umi(data,group,stat):  
    if stat=='uniq':
        umi_counts = data.groupby(group)['umi'].nunique().reset_index()  

    elif stat=='sum':
        umi_counts = data.groupby(group)['umi'].sum().reset_index() 
    group.append('umi')
    umi_counts.columns = group
    return umi_counts

def map_cell(df,meta):
    mapping_dict = meta.set_index('spbc_seq').T.to_dict('list')   
    df[['cell', 'cell_x', 'cell_y']] = df['spbc_seq'].apply(lambda x: pd.Series(mapping_dict.get(x, [None, None, None])))
    #cr_df.shape
    df = df.dropna()
    #df.shape
    return df

def clonebc_umi_percent_percell(df,thr_low=0.1,thr_high=0.9,umi_thr_low=2,umi_thr_high=10
                                   ):
    get_hist_threshold_percent_plot(df,'umi',bins=10,threshold_low=umi_thr_low,threshold_high=umi_thr_high,addlog=True,save=False)
    df1=df[(df['umi']>umi_thr_low) & (df['umi']<umi_thr_high)]
    #print(df1['clonebc'].nunique())
    df1['umi_percent_cell']=df1.groupby(['cell','cell_x','cell_y'])['umi'].transform(lambda x: x / x.sum())#agg(umi_passread_percent=('umi_passed_nread', lambda x: x.sum() / len(x))).reset_index() 
    get_hist_threshold_percent_plot(df1,'umi_percent_cell',bins=10,threshold_low=thr_low,threshold_high=thr_high,save=False)
    #if the clonebc's confident umi percent < 50%, we can consider the clonebc in this cell is contamination
    df1['confident_clonebc']=((df1['umi_percent_cell']>thr_low) & (df1['umi_percent_cell']<thr_high))
    
    return df1
    
def cell_conf_clonebc_percent(df,thr=0.5):
    raw_df=df
    df=df[['cell','cell_x','cell_y','confident_clonebc']].drop_duplicates()
    df=df.groupby(['cell','cell_x','cell_y']).agg(conf_clonebc_percent=('confident_clonebc', lambda x: x.sum() / len(x))).reset_index() 
    get_hist_threshold_percent_plot(df,'conf_clonebc_percent',bins=100,threshold_low=thr,save=False)
    #if the clonebc's confident umi percent < 50%, we can consider the clonebc in this cell is contamination
    df['confident_cell']=df['conf_clonebc_percent']>thr
    get_scatter_plot(df,x='cell_x',y='cell_y',group_col='conf_clonebc_percent',saveplot=False)
    raw_df['confident_cell']=df.set_index('cell')['confident_cell'].reindex(raw_df['cell']).values
    return raw_df

def clonebc_conf_cell_percent(df,thr=0.5):
    raw_df=df
    df=df[['clonebc','cell','confident_cell']].drop_duplicates()
    df= df.groupby('clonebc').agg(clonebc_conf_cell_percent=('confident_cell', lambda x: x.sum() / len(x))).reset_index() 

    get_hist_threshold_percent_plot(df,'clonebc_conf_cell_percent',bins=100,addlog=False,threshold_low=thr,save=False)
    #if the clonebc's confident umi percent < 50%, we can consider the clonebc in this cell is contamination
    df['confident_clonebc_by_cell']=df['clonebc_conf_cell_percent']>thr
    raw_df['confident_clonebc_by_truecell_percent']=df.set_index('clonebc')['confident_clonebc_by_cell'].reindex(raw_df['clonebc']).values
        
    return raw_df