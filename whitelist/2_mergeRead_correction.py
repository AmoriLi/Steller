#2603
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from venn import venn
import random
def hamming_distance(s1, s2):  
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))  

# Function to prepare Lorenz curve data
def lorenz_prep(df, column):
    # Sort the DataFrame by the specified column
    sorted_df = df.sort_values(by=column)
    # Calculate cumulative proportions
    sorted_df['cumulative_population'] = np.linspace(1 / len(sorted_df), 1, len(sorted_df))
    sorted_df['cumulative_income'] = sorted_df[column].cumsum() / sorted_df[column].sum()
    return sorted_df

# Function to calculate Gini index
def gini(df):
    n = len(df)
    B = np.trapz(df['cumulative_income'], df['cumulative_population'])
    A = 0.5 - B
    return A / 0.5

wd=sys.argv[1] #'~/NPC_project/P045_whitelist/w2/'

df1=pd.read_csv(os.path.join(wd,"AA_read.txt"),sep='\t',header=None)
df1=df1.iloc[:, 0].value_counts().reset_index()#df1.groupby(df1.columns[0]).size().reset_index(name='read')
df1.columns=['seq','read']
df1.head()
plt.hist(np.log2(df1['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)')  
plt.ylabel('Frequency') 
plt.savefig(os.path.join(wd,"AA_CloneBC_read_hist_before_collapse.png"))
plt.show()

df2=pd.read_csv(os.path.join(wd,"AB_read.txt"),sep='\t',header=None)
df2=df2.iloc[:, 0].value_counts().reset_index()#df1.groupby(df1.columns[0]).size().reset_index(name='read')
df2.columns=['seq','read']
df2.head()

plt.hist(np.log2(df2['read']), bins=10, edgecolor='black') 
plt.title('Histogram of read')  
plt.xlabel('log2(read)') 
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"AB_CloneBC_read_hist_before_collapse.png"))
plt.show()

df3=pd.read_csv(os.path.join(wd,"BB_read.txt"),sep='\t',header=None)
df3=df3.iloc[:, 0].value_counts().reset_index()#df1.groupby(df1.columns[0]).size().reset_index(name='read')
df3.columns=['seq','read']
df3.head()

plt.hist(np.log2(df3['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)') 
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"BB_CloneBC_read_hist_before_collapse.png"))
plt.show()

df4=pd.read_csv(os.path.join(wd,"BA_read.txt"),sep='\t',header=None)
df4=df4.iloc[:, 0].value_counts().reset_index()#df1.groupby(df1.columns[0]).size().reset_index(name='read')
df4.columns=['seq','read']
df4.head()

plt.hist(np.log2(df4['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)') 
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"BA_CloneBC_read_hist_before_collapse.png"))
plt.show()

df5=pd.read_csv(os.path.join(wd,"trucated_read.txt"),sep='\t',header=None)
df5=df5.iloc[:, 0].value_counts().reset_index()#df1.groupby(df1.columns[0]).size().reset_index(name='read')
df5.columns=['seq','read']
df5.head()

plt.hist(np.log2(df5['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read') 
plt.xlabel('log2(read)') 
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"trucated_CloneBC_read_hist_before_collapse.png"))
plt.show()

df1_set = set(df1['seq'].dropna())
df2_set = set(df2['seq'].dropna())
df3_set = set(df3['seq'].dropna())
df4_set = set(df4['seq'].dropna())

#check the overlap ratio for assessment of library sampling and sequencing saturation
data_dict = {
    "AA": df1_set,
    "AB": df2_set,
    "BB": df3_set,
    "BA": df4_set
}

venn(data_dict)
plt.title("Venn Diagram for 4 sampling")
plt.show()

#A good library whilist sequence read will show obvious normal distribution at log2 scale
df=pd.concat([df1,df2,df3,df4,df5],axis=0).groupby('seq')['read'].sum().reset_index()
df.head
plt.hist(np.log2(df['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)')  
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"Total_CloneBC_read_hist_before_collapse.png"))
plt.show()

#Optional: read filtering based on histogram to remove sequence and PCR noise
df_s=df[df['read']>1]
plt.hist(np.log2(df_s['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)')  
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"filtered_CloneBC_read_hist_before_collapse.png"))
plt.show()

df1_set2 = set(df1['seq'].dropna()) & set(df_s['seq'].dropna())
df2_set2 = set(df2['seq'].dropna()) & set(df_s['seq'].dropna())
df3_set2 = set(df3['seq'].dropna()) & set(df_s['seq'].dropna())
df4_set2 = set(df4['seq'].dropna()) & set(df_s['seq'].dropna())

data_dict2 = {
    "AA": df1_set2,
    "AB": df2_set2,
    "BB": df3_set2,
    "BA": df4_set2
}

venn(data_dict2)
plt.title("Venn Diagram for 4 sampling")
plt.show()

df_s.to_csv(os.path.join(wd,"read_filtered.clonebc.txt"),sep='\t',index=None,header=None)

### Run CloneBC_correction_byStarcode.sh in shell. Afte done, follow next steps:
df2=pd.read_csv(os.path.join(wd,"read_filtered.clonebc.collapsed.txt"),sep='\t',header=None)
df2=df2[[0,1]]
df2.columns=['seq','read']
df2[:5]

plt.hist(np.log2(df2['read']), bins=10, edgecolor='black')  
plt.title('Histogram of read')  
plt.xlabel('log2(read)')  
plt.ylabel('Frequency')  
plt.savefig(os.path.join(wd,"Total_CloneBC_read_hist_after_collapse.png"))
plt.show()

df_set = set(df2['seq'].dropna())

s1 = df_set.intersection(df1_set)
s2 = df_set.intersection(df2_set)
s3 = df_set.intersection(df3_set)
s4 = df_set.intersection(df4_set)
#
data_dict = {
    "AA": s1,
    "AB": s2,
    "BB": s3,
    "BA": s4
}

venn(data_dict)
plt.title("Venn Diagram for 4 sampling")
plt.show()

### HD of random CloneBC sequences
hd = []
for n in range(0,10):
    sampled_bc = list(df2['seq'].sample(n=1000,replace=False).values)
    hamming_distances = []
    for i in range(len(sampled_bc)):
        for j in range(i + 1, len(sampled_bc)):
            # 
            distance = hamming_distance(sampled_bc[i], sampled_bc[j])
            hamming_distances.append(distance)
    if hd is None:
        hd = []
    hd += hamming_distances

plt.hist(hd, bins=(np.max(hd)-np.min(hd)), edgecolor='black', alpha=0.7)
plt.title('Distribution of Hamming Distances')
plt.xlabel('Hamming Distance')
plt.ylabel('Frequency')
plt.grid(False)
plt.savefig(os.path.join(wd,"s1k.HD.png"))
plt.show()

# Prepare Lorenz curve data
df_lorenz = lorenz_prep(df2, 'read')

# Calculate Gini index
gini_index = gini(df_lorenz)
gini_index = f"Gini Index = {gini_index:.2f}"

# Plot Lorenz curve
plt.figure(figsize=(6, 6))
plt.plot(df_lorenz['cumulative_population'], df_lorenz['cumulative_income'], label='Lorenz Curve')
plt.plot([0, 1], [0, 1], 'k--', label='Line of Perfect Equality')
plt.xlabel('Cumulative proportion of CloneBC')
plt.ylabel('Cumulative read of CloneBC')
plt.title('Lorenz Curve and Gini Index')
plt.text(0.1,0.65,gini_index,fontsize=15)
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(wd,"lorenze_curve.png"))
plt.show()

len(np.unique(df2['seq']))

df2.to_csv(os.path.join(wd,"whitelist.txt"),sep="\t",header=None,index=None)

###remove high proportion clonebc
df_merge =df2.copy()
m=np.mean(df_merge['read'])
d=np.std(df_merge['read'])
high_thr=m+4*d
high_thr

#plot the read distribution
#df=df1.copy()
df_merge['percentage'] = df_merge['read'] / df_merge['read'].sum() * 100
#df['group']=df['read'].apply(lambda x: 'high' if x > high_thr else 'normal')
tmp=df_merge.sample(n=10000,replace=False)
tmp=tmp.sort_values("percentage").reset_index()

fig, ax = plt.subplots()
colors = ['red' if value > high_thr else 'grey' for value in tmp['read']]
bars = ax.bar(tmp.index, tmp['read'], color=colors)
ax.set_title('CloneBC read distribution')
ax.set_xlabel('CloneBC')
ax.set_ylabel('Percentage')
plt.show()

# Prepare Lorenz curve data after read cutoff
df_lorenz = lorenz_prep(df_merge[df_merge['read']<=high_thr], 'read')

# Calculate Gini index
gini_index = gini(df_lorenz)
gini_index = f"Gini Index = {gini_index:.2f}"

plt.figure(figsize=(6, 6))
plt.plot(df_lorenz['cumulative_population'], df_lorenz['cumulative_income'], label='Lorenz Curve')
plt.plot([0, 1], [0, 1], 'k--', label='Line of Perfect Equality')
plt.xlabel('Cumulative proportion of CloneBC')
plt.ylabel('Cumulative read of CloneBC')
plt.title('Lorenz Curve and Gini Index')
plt.text(0.1,0.65,gini_index,fontsize=15)
plt.legend()
plt.grid(True)
#plt.savefig(os.path.join(wd,"lorenze_curve.png"))
plt.show()

df_filt=df_merge[df_merge['read']<=high_thr]
df_filt['seq'].nunique()
df_filt.to_csv(os.path.join(wd,"CloneBC.withoutHighProp.txt"),sep="\t",header=None,index=None)
