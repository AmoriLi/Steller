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
# Use the Jaccard similarity formula based on intersections and unions  
def jaccard_similarity(X):  
    # Calculate the intersection, which is the sum of min of the rows  
    intersection = X.dot(X.T).A  
    # Calculate the union  
    row_sums = X.sum(axis=1).A.ravel()  
    union = row_sums[:, None] + row_sums[None, :] - intersection  
    # Return Jaccard similarity  
    return intersection / union  

# Function to call clones  
def call_clones(jac_mat_low, cells, jac_th=0.6, return_graph=False):  
    '''  
    Performs cellBC clone calling and returns a clone table.  
    - jac_mat_low: Jaccard similarity matrix  
    - cells: list or array of cell barcodes  
    - jac_th: Jaccard threshold for clone calling  
    - return_graph: flag to return the graph object  
    '''  
    # Create a binary adjacency matrix based on the Jaccard threshold  
    jac_mat = jac_mat_low > jac_th  

    # Create edges DataFrame from non-zero entries in the adjacency matrix  
    edge_df = pd.DataFrame(np.array(jac_mat.nonzero()).T, columns=['source', 'target'])  
    
    # Preparing the vertex DataFrame  
    vertex_df = pd.DataFrame({'idx': range(len(cells)), 'cb': cells})  
    vertex_df['cb'] = vertex_df['cb'].astype(str)  # Ensure cell barcodes are treated as strings  

    # Create the clone graph using igraph  
    g = ig.Graph.DataFrame(edges=edge_df, directed=False, use_vids=False, vertices=vertex_df)  
    
    # Remove isolated vertices  
    g.vs.select(_degree=0).delete()  

    # Call clones using connected components  
    clones = g.components()  

    # Evaluate columns for clone table  
    edge_den = [g.induced_subgraph(i).density() for i in clones]  
    clones_bc = [g.induced_subgraph(i).vs['cb'] for i in clones]  
    clone_id = [[j + 1] * len(i) for j, i in enumerate(clones_bc)]  
    
    # Create the clone table as a DataFrame  
    clone_table = pd.concat([pd.DataFrame({'clone.id': i, 'cell': j, 'edge.den': k})   
                             for i, j, k in zip(clone_id, clones_bc, edge_den)], ignore_index=True)  

    # Filter out clones that have only one cell (single-cell clones)  
    clone_table = clone_table.groupby('clone.id').filter(lambda x: len(x) > 1)  

    if return_graph:  
        return g, clone_table  

    return clone_table  # if return_graph is False  

### 2026-02: debug max(edge_den)>1
def call_clones_v2(jac_mat_low, cells, jac_th=0.6, return_graph=False):
    jac_mat = jac_mat_low > jac_th
    
    jac_mat_triu = np.triu(jac_mat, k=1)

    edges = np.array(jac_mat_triu.nonzero()).T
    edge_df = pd.DataFrame(edges, columns=['source', 'target'])
    vertex_df = pd.DataFrame({'cb': cells})
    vertex_df['idx'] = vertex_df.index 

    vertex_df = vertex_df[['idx', 'cb']]
    
    g = ig.Graph.DataFrame(edges=edge_df, directed=False, vertices=vertex_df)

    g.vs.select(_degree=0).delete()

    clones = g.connected_components()
    all_rows = []
    for j, nodes in enumerate(clones):
        subgraph = g.induced_subgraph(nodes)
        density = subgraph.density()
        cb_list = subgraph.vs['cb']
        
        sub_df = pd.DataFrame({
            'clone.id': j + 1,
            'cell': cb_list,
            'edge.den': density
        })
        all_rows.append(sub_df)
    
    if not all_rows:
        return (g, pd.DataFrame()) if return_graph else pd.DataFrame()
        
    clone_table = pd.concat(all_rows, ignore_index=True)

    if return_graph:
        return g, clone_table
        
    return clone_table

def plot_jaccard_heatmap(jaccard_matrix, cells, clone_table,wd):  
    # Filter clone_table to only include relevant cells  
    relevant_cells = clone_table['cell'].values  

    # Count the number of cells in each clone  
    clone_sizes = clone_table['clone.id'].value_counts()  

    # Create a DataFrame that maps each cell to its corresponding clone size  
    size_map = clone_table.groupby('cell')['clone.id'].first().map(clone_sizes)  

    # Sort relevant cells by their clone sizes, and get the sorted indices  
    sorted_cell_names = size_map.sort_values(ascending=False).index.values 
    #sorted_cell_names = relevant_cells[sorted_cells]  

    # Generate filtered and sorted Jaccard matrix  
    relevant_indices = [np.where(cells == cell)[0][0] for cell in sorted_cell_names if cell in cells]  
    filtered_jaccard_matrix = jaccard_matrix[relevant_indices, :][:, relevant_indices]  
    # Set the diagonal to 1  
    np.fill_diagonal(filtered_jaccard_matrix, 1) 
    # Create a new DataFrame for the heatmap  
    heatmap_df = pd.DataFrame(filtered_jaccard_matrix, index=sorted_cell_names, columns=sorted_cell_names)  

    # Create the heatmap  
    plt.figure(figsize=(10, 8))  
    sns.heatmap(heatmap_df, vmin=0,vmax=1,xticklabels=False,yticklabels=False, annot=False, cmap="Blues", square=True, linewidths=0)  
    plt.title('Jaccard Similarity Heatmap of Cells with Clones')  
    plt.xlabel('Cells')  
    plt.ylabel('Cells')  
    #plt.xticks(rotation=45)  
    #plt.yticks(rotation=45)  
    plt.tight_layout()  
    plt.show()  
    plt.savefig(os.path.join(wd,"clone.jc.heatmap.pdf"))

def intBC_heatmap_modified(allele_table, group ='clone.id',sample_n=300,wd):
    unique_vals = allele_table.loc[:,group].unique()  # 获取唯一值
    sampled_vals = random.sample(list(unique_vals), min(sample_n, len(unique_vals)))
    sampled_df = allele_table[allele_table.loc[:,group].isin(sampled_vals)].copy()
    allele_table=sampled_df
    at_pivot_I = pd.pivot_table(
                allele_table,
                index="cell",
                columns="clonebc",
                values="umi",
                aggfunc="sum",
                fill_value=0  # 避免NaN
            )
    plt.close()
    flat_master = []
    flat_cellBC = []
    es_annotation = []
    for n, lg in allele_table.groupby(group):
        for item in lg["clonebc"].unique():
            flat_master.append(item)
        for item in lg["cell"].unique():
            flat_cellBC.append(item)
        es_annotation += [n]*len(lg['cell'].unique())
    
    new_flat_master =  list({}.fromkeys(flat_master).keys()) # remove duplicate intBC

    at_pivot_I = at_pivot_I.loc[flat_cellBC,new_flat_master]
    
    colors = list(plt.get_cmap('tab10').colors)
    module_colors = {j: colors[(i-1) % len(colors)] for i, j in zip(range(len(np.unique(es_annotation))),
                                                                    np.unique(es_annotation))}

    row_colors = pd.Series([module_colors[i] for i in es_annotation], index = at_pivot_I.index)
    cbar_mat = np.zeros((len(row_colors), 1, 3))
    for r in range(len(row_colors)):
        cbar_mat[r, 0, :] = row_colors.iloc[r][0], row_colors.iloc[r][1], row_colors.iloc[r][2]
    h2 = plt.figure(2,figsize=(7,7))
    axmat2 = h2.add_axes([0.3,0.1,0.6,0.8])
    cmap = matplotlib.colormaps['gray_r']
    im2 = axmat2.matshow(at_pivot_I, cmap = cmap, vmin=0, vmax=10, aspect='auto', origin='upper')
    # im2 = axmat2.matshow(at_pivot_I, aspect='auto', origin='upper')
    axmat2.set_yticks([])
    axmat2.set_xticks([])
    cbar = h2.add_axes([0.2, 0.1, 0.05, 0.8])
    im3 = cbar.imshow(cbar_mat, aspect='auto', origin='upper', interpolation="none")
    plt.title('Consensus intBC pivot table')
    plt.savefig(os.path.join(wd,'clonebc_pivot_heatmap.png'))