# scRNA-seq-breast-cancer

In this project, we have analysed single cell RNA sequencing data collected from human breast cancer samples from a published study GSE225600 (Liu et al., 2023). Breast cancer is a complex and heterogeneous disease comprising of varied cellular subtypes that play a role in the manifestation of the disease. Every cell undergoes differential regulation mechanisms to express different cellular phenotypes such as cell markers or proteins. single-cell RNA sequencing has emerged as a powerful tool to understand this cellular heterogeneity. In this report, the computational analysis of single-cell omics between primary tumors and metastatic tumors of breast cancer patients has been performed. We visualise the apparent cellular heterogeneity of the tissue samples through state-of-the-art dimensional reduction techniques such as t-SNE and UMAP. We further map out the clusters using clustering algorithms on the reduced dimensional representation and annotate these clusters through reference datasets. We further study the cellular diversity through state-of-the-art statistical methods using over-representation analysis, gene set enrichment analysis and compositional studies. We found out that there were important cell types that varied in both conditions such as endothelial cells and epithelial cells. We also found epithelial-mesenchymal transition and N-F- Kappa B signalling pathway to be differentially expressed in both the conditions. Overall, the computational analysis of single-cell RNA sequencing has been extensively discussed and used to understand the molecular and cellular complexity of breast cancer. This technology may contribute to an increased understanding of the intricate mechanisms that take place for a cancer to be more invasive. This project was undertaken as part of the Indian Academy of Sciences summer fellowship programme.

Dependencies:

1. preprocessing.py : It's used to preprocess the raw files generated using Cell Ranger's count function. This also generates QC plots for the cells and the genes in the dataset.
2. soupXcounts_per_sample.ipynb : It's a Jupyter notebook that has various functionalities such as sample preprocessing, normalization, dimensionality reduction, SoupX contamination rate, doublet finder and clustering.
3. soupXcounts_all_samples.ipynb : It's a Jupyter notebook that concatenates different samples and shows QC metrics across the various samples. In addition, it shows the dimensionality reduction plots coloured by samples to show the distribution of samples.
4. sample-integration.py : This file has several batch correction methods that can be used for generating corrected .h5ad files.
5. annotation.ipynb : It's a Jupyter notebook that annotates the clusters through a logistic-regression based method called CellTypist.
6. enrichment_compositional.ipynb : This file has all the differential expression analysis related codes where both overexpression analysis and functional scoring for individual pathways have been performed.
7. compositional_cell.py : This file uses SCCODA model for performing compositional analysis on cell types.

Current system info generated using R:
   		 R version 4.3.1 (2023-06-16)
   		 Platform: x86_64-pc-linux-gnu (64-bit)
   		 Running under: Ubuntu 20.04.6 LTS

   		 Matrix products: default
   		 BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
   		 LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

   		 locale:
   		  [1] LC_CTYPE=en_IN   	LC_NUMERIC=C     	LC_TIME=en_IN  	 
   		  [4] LC_COLLATE=en_IN 	LC_MONETARY=en_IN	LC_MESSAGES=en_IN   
   		  [7] LC_PAPER=en_IN   	LC_NAME=C        	LC_ADDRESS=C   	 
   		 [10] LC_TELEPHONE=C   	LC_MEASUREMENT=en_IN LC_IDENTIFICATION=C

   		 time zone: Asia/Kolkata
   		 tzcode source: system (glibc)

   		 attached base packages:
   		 [1] stats4	tools 	stats 	graphics  grDevices utils 	datasets
   		 [8] methods   base	 

   		 other attached packages:
   		  [1] scry_1.12.0             	BiocParallel_1.34.2   	 
   		  [3] scDblFinder_1.14.0      	scater_1.28.0         	 
   		  [5] ggplot2_3.4.2           	scuttle_1.10.1        	 
   		  [7] SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2
   		  [9] Biobase_2.60.0          	GenomicRanges_1.52.0  	 
   		 [11] GenomeInfoDb_1.36.1     	IRanges_2.34.1        	 
   		 [13] S4Vectors_0.38.1        	BiocGenerics_0.46.0   	 
   		 [15] MatrixGenerics_1.12.2   	matrixStats_1.0.0     	 
   		 [17] SeuratObject_4.1.3      	Seurat_4.3.0.1        	 
   		 [19] Matrix_1.5-4.1          	SoupX_1.6.2           

Steps that were followed to install the dependencies:

1. First of all, use the requirements.txt to install the dependencies using:
   	 pip install -r requirements.txt
This will install the python based requirements.

2. To use R packages with Python's data objects, a package called rpy2 is required. This pipeline used rpy2 (version 3.5.8) 

3. While installing rpy2, make sure you already have the latest version of R on your pc. At the moment, this pipeline was written using R (version 4.3.1)

4. To show the plots in VSCode, matplotlib requires tkinter so make sure it is downloaded on the Ubuntu platform.

5. R packages are installed separately either using RStudio or in the Jupyter notebook using %%R. At times, extra installations using apt-install have been performed for installing libxml2-dev, xml2-config, libtiff5-dev and libxt-dev which were dependencies for scDblFinder and scrat.

6. BBKNN, pertpy and leiden algorithm also needs to be installed separately for sample_integration.py and cell composition analysis. These requirements are present in the requirements.txt file.


