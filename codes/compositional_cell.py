import pandas as pd

import scanpy as sc
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pertpy as pt # used for sccoda model

print("Running...")
adata_1 = sc.read_h5ad("/home/digemed/Sunaina-single-cell/samples/lymph_node_nbiswbatch_corrected.h5ad")
adata_2 = sc.read_h5ad("/home/digemed/Sunaina-single-cell/samples/primary_tumor_nbiswbatch_corrected.h5ad")

adata_1.obs['condition'] = 'Metastatic'
adata_2.obs['condition'] = 'Primary' #setting condition for both the adata

adata = adata_1.concatenate(adata_2)

print(adata.obs['sample'].value_counts())
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    adata, type="cell_level", generate_sample_level=True, cell_type_identifier="majority_voting", sample_identifier="sample",
    covariate_obs=['condition']
) # loading the data into the sccoda model which prints the covariate and intercept values for different cell types
print(sccoda_data)

pt.pl.coda.stacked_barplot(
    sccoda_data, modality_key="coda", feature_name="condition",
)
plt.tight_layout()
plt.savefig("stacked_bar.png", bbox_inches='tight') #stacked bar plot showing the compositional changes

sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda",
    formula="condition",
    reference_cell_type="automatic",
)

sccoda_model.run_nuts(sccoda_data, modality_key="coda", rng_key=1234)

print(sccoda_data["coda"].varm['effect_df_condition[T.Primary]'].to_string())
print(sccoda_data["coda"].varm['intercept_df'].to_string())

sccoda_model.set_fdr(sccoda_data, 0.05) #setting FDR to be 0.05 to be less stringent
sccoda_model.credible_effects(sccoda_data, modality_key="coda")

print(sccoda_data["coda"].varm['effect_df_condition[T.Primary]'].to_string())
print(sccoda_data["coda"].varm['intercept_df'].to_string())

pt.pl.coda.effects_barplot(sccoda_data, "coda", "condition")
plt.savefig('fold_changes_fdr_0.05.png', bbox_inches='tight') #final foldchange values according to the sccoda model

sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)

sc.pl.umap(adata, color=["condition", "sample", "majority_voting"], ncols=3, wspace=1) #showing UMAP embeddings coloured by different parameters for all the samples
plt.savefig('UMAP.png', bbox_inches='tight')
