import scanpy as sc
import pandas as pd
import sccoda


adata_1 = sc.read_h5ad("/home/digemed/Sunaina-single-cell/samples/lymph_node_nbiswbatch_corrected.h5ad")
adata_2 = sc.read_h5ad("/home/digemed/Sunaina-single-cell/samples/primary_tumor_nbiswbatch_corrected.h5ad")

adata_1.obs['condition'] = 'Metastatic'
adata_2.obs['condition'] = 'Primary'

adata = adata_1.concatenate(adata_2)
df = adata.to_df()
df['majority_voting'] = adata.obs['majority_voting']
cell_types = df.groupby('majority_voting').mean()

