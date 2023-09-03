import pandas as pd
import scanpy as sc
import numpy as np

def initiation(mtx_file, features):
    """ 
    This function takes in the raw counts matrix and initialises it as an AnnData data object. Further, it also takes in annotations that are used to 
    provide index to both the observation and variable rows.

    Parameters: 
        1. mtx_file: Counts matrix in the matrix format  
        2. feature_dir: features directory

    """

    adata = sc.read_10x_mtx(mtx_file, cache= True) #reading the matrix file (.mtx) as an annData object
    features.columns = ['ENSEMBL', 'GENESYM', 'TYPE'] # providing column names to the features
    # adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)] # Providing IDs to the cells: Cell1, Cell2, ... and so on.
    adata.var_names = list(features['GENESYM']) # Prodiving gene names as IDs
    adata.var_names_make_unique() # say ABC12 is a duplicate, then it'd simply give ABC12.1, ABC12.2 as the label
    return adata

def quality_check(adata, n_genes_counts = 10000, pct_mt = 10, u_min_genes = 200):
    """ 
    This function takes in the annData data object and calculates the amount of mitochondrial, ribosome and haemoglobin genes. 
    n_genes_counts were filtered according to: (https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)

    Parameters: 
        1. adata
        2. n_genes_counts is the threshold of read counts. Extremely high number of counts are associated with doublets. Extremely low number of counts are associated with dying cells.
        3. pct_mt is the threshold for percentage of mitochondrial genes that a particular cell has. 

    """

    adata.var['mt'] = adata.var_names.str.startswith("MT-") # mitochondrial genes have TRUE in this column, others have FALSE
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS", "RPL")) #ribosomal genes
    adata.var["haemoglobin"] = adata.var_names.str.contains(("^HB[^(P)]")) #haemoglobin genes
    adata.var_names_make_unique() # it makes the variable names unique: ABC1.1, ABC1.2 and so on.
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "haemoglobin"], inplace=True, percent_top=[20], log1p=True) # calculates QC metrics such as percentage counts in top 20 genes or total counts of haemoglobin genes
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', "pct_counts_haemoglobin", "pct_counts_ribo"], multi_panel=True) # any QC metrics can be explored here as the violin plot

    adata = adata[adata.obs.n_genes_by_counts < n_genes_counts, :] # if you want to keep the gene count below  a range
    adata = adata[adata.obs.pct_counts_mt < pct_mt, :] # 

    sc.pp.filter_cells(adata, min_genes= u_min_genes) # removing cells that don't have minimum 200         
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color = "pct_counts_mt") 
    return adata

def cleaning_the_soup(adata, raw_dir):

    """
    Correcting for ambient mRNA. This method requires the clusters of the cells to be already initialized. Thus, the annData object is first normalized, its PCA embedding calculated so that
    the lower dimension neighbor embedding is formed. Then, Leiden clustering was performed. This embedding is further used to correct the contamination fraction.

    Parameters: 
    1. adata
    2. raw_dir = directory for the raw matrix counts

    """

    adata_pp = adata.copy() 
    adata_pp.var_names_make_unique()
    adata.var_names_make_unique()

    sc.pp.normalize_per_cell(adata_pp) 
    sc.pp.log1p(adata_pp) # log shift normalization explained later
    sc.pp.pca(adata_pp) # PCA embedding

    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added= "soupx_groups" )

    # making the parameters ready for soupx groups
    soupx_groups = adata_pp.obs["soupx_groups"]
    cells = adata.obs_names
    genes = adata.var_names
    
    data = adata.X.T

    raw_matrix = sc.read_10x_mtx(raw_dir, cache=True)
    raw_matrix.var_names_make_unique()
    data_raw = raw_matrix.X.T
    
    return data, data_raw, genes, cells, soupx_groups # filtered matrix, raw matrix, list of genes, list of cells and leiden clusters named as soupx_groups
    

