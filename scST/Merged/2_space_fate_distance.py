import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from matplotlib.lines import Line2D
wd='~/NPC_project/BMK/tissue/brain/B33_4_P14_half/section1/Merged_2512'
outdir=os.path.join(wd,'4_lineageCP')
df = pd.read_csv(os.path.join(wd,'3_clone/test7a.snST.scST.integrated.final_ct23.L123_region.umap.csv'),index_col=0)
df[:5]
df1 = df[df['clone.id2'] != 'none'][df['clone.id2'].notna()].copy()
mouse = 'A13'
df1 = df1[df1['mouse']==mouse][~df1['final_ct2'].isin(['Mixed','Other','Myeloid','Unknown'])]
clone_counts = df1['clone.id2'].value_counts()
selected_clones = clone_counts[clone_counts >= 2].index
df2 = df1[df1['clone.id2'].isin(selected_clones)].copy()

coords_xy = df2[['x_ad', 'y_ad']].values
coords_umap = df2[['UMAP_1', 'UMAP_2']].values
clones = df2['clone.id2'].values

dist_xy_all = pdist(coords_xy)
dist_umap_all = pdist(coords_umap)

dist_xy_within = []
dist_umap_within = []

for clone in np.unique(clones):
    mask = (clones == clone)
    if mask.sum() > 1:
        dist_xy_within.extend(pdist(coords_xy[mask]))
        dist_umap_within.extend(pdist(coords_umap[mask]))

dist_xy_within = np.array(dist_xy_within)
dist_umap_within = np.array(dist_umap_within)

plt.figure(figsize=(6,3))

# spatial distance
plt.subplot(1, 2, 1)
plt.ecdf(dist_xy_within, label='Within a Clone')
plt.ecdf(dist_xy_all, label='All Cells (Background)')
plt.title('Spatial Distance (XY)')
plt.xlabel('Distance (um)')
plt.ylabel('ECDF')
plt.legend()

#  UMAP fate distance
plt.subplot(1, 2, 2)
plt.ecdf(dist_umap_within, label='Within a Clone')
plt.ecdf(dist_umap_all, label='All Cells (Background)')
plt.title('Transcriptomic Distance (UMAP)')
plt.xlabel('UMAP distance')
plt.ylabel('ECDF')
plt.legend()

plt.tight_layout()
plt.savefig(os.path.join(outdir,mouse+'.spatial.umap.distance.within_clone.pdf'))
plt.show()

### spatial distance within major fate direction
df2=pd.read_csv(os.path.join(wd,'3_clone/test7a.scST.clone_cell.lineage_final_ct2.error_adjusted.meta.csv'),index_col=0)
df2=df2[df2['mouse']=='A13']
df2=df2[df2['fate']!='Undefined']
df2[:5]
clone_counts = df2['clone'].value_counts()
selected_clones = clone_counts[clone_counts >= 2].index
df3 = df2[df2['clone'].isin(selected_clones)].copy()
results = []

for cell_type in df3['fate'].unique():
    subset = df3[df3['fate'] == cell_type]
    type_clones = subset['clone'].unique()
    
    for clone in type_clones:
        clone_subset = subset[subset['clone'] == clone]
        
        if len(clone_subset) >= 2:
           
            coords = clone_subset[['x_ad', 'y_ad']].values
            distances = pdist(coords)
            
            for d in distances:
                results.append({
                    'CellType': cell_type,
                    'Distance': d
                })
                
dist_df = pd.DataFrame(results)

min_observations = 2
counts = dist_df['CellType'].value_counts()
keep_types = counts[counts >= min_observations].index
dist_df = dist_df[dist_df['CellType'].isin(keep_types)].copy()

num_random_samples = min(1000, len(df3)) 
random_cells_coords = df2[['x_ad', 'y_ad']].sample(n=num_random_samples, random_state=42).values
background_distances = pdist(random_cells_coords)

stats = dist_df.groupby('CellType')['Distance'].median().sort_values()
threshold = 2000

clustered_types = stats[stats <= threshold].index.tolist()
dispersed_types = stats[stats > threshold].index.tolist()

colors_c = sns.color_palette("Reds_r", len(clustered_types)).as_hex()
colors_d = sns.color_palette("Blues_r", len(dispersed_types)).as_hex()

palette_dict = {t: colors_c[i] for i, t in enumerate(clustered_types)}
palette_dict.update({t: colors_d[i] for i, t in enumerate(dispersed_types)})

plt.figure(figsize=(6, 4))
ax = plt.gca()

sns.ecdfplot(background_distances, label='Random Background', color='grey', linestyle='--', linewidth=2.5, ax=ax, alpha=0.8)

for c_type in stats.index:
    data_subset = dist_df[dist_df['CellType'] == c_type]['Distance']
    sns.ecdfplot(data_subset, label=c_type, color=palette_dict[c_type], linewidth=2, ax=ax, alpha=0.8)

legend_elements = [Line2D([0], [0], color='grey', linestyle='--', lw=2, label='Random Background')]
legend_elements.append(Line2D([0], [0], color='none', label=' ')) 

legend_elements.append(Line2D([0], [0], color='none', label=f'$\mathbf{{NEAR\ (\leq {threshold}um)}}$'))
for t in clustered_types:
    legend_elements.append(Line2D([0], [0], color=palette_dict[t], lw=2, label=t))

legend_elements.append(Line2D([0], [0], color='none', label=' ')) 
legend_elements.append(Line2D([0], [0], color='none', label=f'$\mathbf{{FAR\ (> {threshold}um)}}$'))
for t in dispersed_types:
    legend_elements.append(Line2D([0], [0], color=palette_dict[t], lw=2, label=t))

ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9, frameon=False)

plt.title(f'Clonal Dispersion with Random Background (Threshold: {threshold} $\mu m$)', fontsize=8, pad=20)
plt.xlabel('Intra-clone Spatial Distance ($\mu m$)')
plt.ylabel('Cumulative Fraction')
plt.grid(True, which="both", ls=":", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(outdir,mouse+'RGC.fate.clonal_dispersion.pdf'))
plt.show()

