import scanpy as sc
import anndata

lymph_3 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525752_GSM7050930_PT_3_lymph_node_sc_S3_normalized.h5ad")
lymph_2 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525752_GSM7050930_PT_2_lymph_node_S2_normalized.h5ad")
lymph_1 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525752_GSM7050930_PT_7_lymph_node_S7_normalized.h5ad")

primary_1 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525754_GSM7050928_PT_7_tumor_sc_S7_normalized.h5ad")
primary_2 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525757_GSM7050925_PT_2_tumor_sc_S2_normalized.h5ad")
primary_3 = anndata.read_h5ad("/home/digemed/Sunaina-single-cell/GSE225600/h5ad samples/SRR23525756_GSM7050926_PT_3_tumor_sc_S3_normalized.h5ad")

lymph_1.var_names_make_unique()
lymph_2.var_names_make_unique()
lymph_3.var_names_make_unique()
primary_1.var_names_make_unique()
primary_2.var_names_make_unique()
primary_3.var_names_make_unique()

lymph_1.obs['sample'] = "Lymph_1"
lymph_2.obs['sample'] = "Lymph_2"
lymph_3.obs['sample'] = "Lymph_3"

primary_1.obs["sample"] = "Primary_1"
primary_2.obs["sample"] = "Primary_2"
primary_3.obs["sample"] = "Primary_3"

primary_adata = primary_1.concatenate(primary_2, primary_3)
adata = primary_adata

adata.var_names_make_unique()
adata.obs_names_make_unique()
print(adata)
print("Total number of observations/cells in each sample: ")
print(adata.obs["sample"].value_counts())
print(adata)
label_key = "majority_voting"
batch_key = "sample"
adata.uns[batch_key + "colors"] = [
    "#5ef160",
    "#0b034b",
    "#c1fe57",
    "#7501bf",
    "#811473"
]


# adata_hvg = adata[:, adata.var["highly_variable"]].copy()


adata_hvg = adata
normalize = 1 #change this parameter accordingly
if normalize == 0: #raw data
    sc.pp.neighbors(adata)
    sc.pp.pca(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[label_key, batch_key], wspace=1)    

elif normalize == 1:
    sc.pp.pca(adata_hvg)
    sc.external.pp.bbknn(adata_hvg, batch_key="sample") #BBKNN method of batch correction
    sc.tl.umap(adata_hvg)
    sc.pl.umap(adata_hvg, color= "sample")


elif normalize == 2:
    sc.pp.combat(adata_hvg, key="sample", inplace= True) #COMBAT method of batch correction
    sc.pp.pca(adata_hvg)
    sc.pp.neighbors(adata_hvg)
    sc.pp.pca(adata_hvg)
    sc.tl.umap(adata_hvg)
    sc.pl.umap(adata_hvg, color=[batch_key], wspace=1)
elif normalize == 3:
    sc.external.pp.mnn_correct(adata_hvg, batch_key= "sample", inplace=True) #MNN method of batch correction
    sc.pp.neighbors(adata_hvg)
    sc.pp.pca(adata_hvg)
    sc.tl.umap(adata_hvg)
    sc.pl.umap(adata_hvg, color=[batch_key], wspace=1)  
else:
    print("")

# 
# sc.pp.highly_variable_genes(
#     adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key
# ) # after performing batch correction, one can find the variable genes while taking the batch information into account.
# adata_hvg.write_h5ad("E:/IAS-SRFP-2023/Single_Cell_Portal/top_2000_batch_corrected.h5ad")
