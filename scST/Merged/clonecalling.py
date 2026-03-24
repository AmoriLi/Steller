import sys
import os
import glob
import pandas as pd
import numpy as np
import scipy
import igraph as ig
from scipy.io import mmread
import matplotlib.colors as mcolors  

from scipy.sparse import csr_matrix, csc_matrix  
from sklearn.metrics import jaccard_score  
from sklearn.metrics import pairwise_distances  

#import adjustText
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.append("~/NPC_project/BMK/pipeline/amplicon")
import LG_utils
sys.setrecursionlimit(10000)

wd=sys.argv[1]#"~/NPC_project/scRNAseq/B12_P6/integrate/2_clonebc"
jc_thr=0.7#sys.argv[2]
outdir=os.path.join(wd,"tmp")

if not os.path.exists(outdir):
    os.makedirs(outdir)

ma=mmread(glob.glob(os.path.join(wd,"*matrix.mtx"))[0]).T 
cells=pd.read_csv(glob.glob(os.path.join(wd,"*cell.meta.tsv"))[0],sep='\t').index.values 
celltags=pd.read_csv(glob.glob(os.path.join(wd,"*clonebc.tsv"))[0],sep='\t',index_col=0,header=None).index.values

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
plt.savefig(os.path.join(wd,"jc.hist.png"))

clones=LG_utils.call_clones_v2(jaccard_matrix,cells,0.7)

clones.to_csv(os.path.join(outdir,'clones.csv'), index=False) 

LG_utils.plot_jaccard_heatmap(jaccard_matrix,cells,clones,wd=wd)

