import sys
import os
import glob
import pandas as pd
import numpy as np
import scipy
import igraph as ig
from scipy.io import mmread
import matplotlib.colors as mcolors  
import matplotlib.patches as mpatches
from scipy.sparse import csr_matrix, csc_matrix  
from sklearn.metrics import jaccard_score  
from sklearn.metrics import pairwise_distances  

#import adjustText
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.append("/public/home/zxli_gibh/NPC_project/BMK/pipeline/amplicon")
import LG_utils
sys.setrecursionlimit(10000)

wd="~/NPC_project/BMK/tissue/brain/P5/Merged_2512/3_clone"
jc_thr=0.7#sys.argv[2]
outdir=os.path.join(wd,"tmp")
if not os.path.exists(outdir):
    os.makedirs(outdir)

ma=mmread(glob.glob(os.path.join(wd,"*matrix.mtx"))[0]).T #need change
cells=pd.read_csv(glob.glob(os.path.join(wd,"*cell.meta.tsv"))[0],sep='\t').index.values #need change
celltags=pd.read_csv(glob.glob(os.path.join(wd,"*clonebc.tsv"))[0],sep='\t',index_col=0,header=None).index.values #need change
#jac_mat = jaccard_similarity(ma.T)
#jac_mat.setdiag(1)
umi_table = pd.read_csv(os.path.join(wd,'cell.clonebc.umi_table.txt'),header=0,index_col=0)

if not isinstance(ma, csr_matrix):  
    ma = ma.tocsr()  

# Convert to binary format for Jaccard similarity  
binary_matrix = (ma > 0).astype(int) 

# Compute the Jaccard similarity matrix  
jaccard_matrix = LG_utils.jaccard_similarity(binary_matrix)  
np.fill_diagonal(jaccard_matrix, 1)
# Calculate some statistics
nonzero_jaccard = jaccard_matrix[np.triu_indices(len(jaccard_matrix), k=1)]
nonzero_jaccard = nonzero_jaccard[nonzero_jaccard > 0]

# Create histogram of non-zero Jaccard indices
plt.figure(figsize=(4, 3))
plt.hist(nonzero_jaccard, bins=50, edgecolor='black')
plt.title('Distribution of Non-zero Jaccard Indices')
plt.xlabel('Jaccard Index')
plt.ylabel('Frequency')
plt.show()

clones = pd.read_csv(os.path.join(outdir,'clones.csv')) 
clones[:5]
clone_sizes = clones['clone.id'].value_counts().sort_values(ascending=False)
clone_sizes
# clone.id by size descending
sorted_clones = clone_sizes.index.values
sorted_clones
sorted_cell_names = []
for clone in sorted_clones:
    cells_in_clone = clones[clones['clone.id'] == clone]['cell'].values
    sorted_cell_names.extend(cells_in_clone)

relevant_indices = [np.where(cells == cell)[0][0] for cell in sorted_cell_names]
    
# filter and sort jaccard matrix
filtered_jaccard_matrix = jaccard_matrix[np.ix_(relevant_indices, relevant_indices)]

#
np.fill_diagonal(filtered_jaccard_matrix, 1)

clone_series = clones.set_index('cell').loc[sorted_cell_names]['clone.id']

#random color for each cloneid
unique_clones = clone_series.unique()
palette = sns.color_palette('tab20', n_colors=len(unique_clones))
clone2color = dict(zip(unique_clones, palette))

g = sns.clustermap(
    pd.DataFrame(filtered_jaccard_matrix, index=sorted_cell_names, columns=sorted_cell_names),
    row_colors=clone_series.map(clone2color),
    #col_colors=lg_series.map(lg2color),
    cmap='Blues', vmin=0, vmax=1,
    linewidths=0, figsize=(10, 10),
    row_cluster=False, col_cluster=False,
    dendrogram_ratio=0,
    xticklabels=False, yticklabels=False,
    cbar_pos=None
)
plt.savefig(os.path.join(wd,'sorted.jac.heatmap.png'),dpi=400,bbox_inches='tight',
           format='png')
plt.show()

### pivot table heatmap
umi_table=pd.read_table(os.path.join(wd,"cell.clonebc.umi_table.txt"),sep=',',index_col=0)
umi_table[:5]
merged_table = pd.merge(umi_table, clones, on='cell', how='inner').drop_duplicates()
merged_table
LG_utils.intBC_heatmap_modified(merged_table,wd=wd)