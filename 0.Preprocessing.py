
import scanpy as sc
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
import rpy2.robjects as robjects
from rpy2.robjects import r
import hdf5plugin
import seaborn as sns
import colorcet as cc
from pathlib import Path
from matplotlib.pyplot import rc_context
from scipy.sparse import csr_matrix, issparse
import matplotlib as mpl
import itertools
import anndata2ri
import scvi
from matplotlib.backends.backend_pdf import PdfPages
import re

robjects.pandas2ri.activate()
anndata2ri.activate()
inferCNV = False
print("si parte!\n LET'S GOSKY")


print("###########################################################")
print("load data")
print("###########################################################")
# file xlsx containing in the first column the sample names and in the other column the clinical/biological metadata (including headers)
# the first column header should always be "orig.ident"
file_path = './ptsinfo.xlsx'

# Leggi il file Excel
df_clinici = pd.read_excel(file_path)

df_clinici.columns
df_clinici = df_clinici.loc[:, ['orig.ident']]
df_clinici = df_clinici.dropna(subset=['orig.ident'])
#directory containing all cellranger count outputs
All_Patients = os.listdir("./cellranger")
to_analysis = All_Patients
Patients = [p for p in All_Patients if p in to_analysis and p != ".DS_Store" and p != "Multiplex"]
Patients.sort()
print(Patients)
### read all files
adata_union = {}
for p in Patients:
    adata = sc.read_10x_mtx(
        "cellranger/"+p+"/outs/filtered_feature_bc_matrix/",  # the directory with the 3 files(barcodes, features, matrix) 
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=True,  # write a cache file for faster subsequent reading
    )
    adata_union[p] = adata

for name, adata in adata_union.items():
    #remove dates from sample names if present
    adata.obs['orig.ident'] = re.sub(r'20\d{2}-', '', name)
    

adata_union_updated = {}

# dictionary containing adata 
for key, adata in adata_union.items():
    #r emove dates from sample names keys if present
    new_key = re.sub(r'20\d{2}-', '', key)
    
    # add the new key and corresponding AnnData to the dictionary
    adata_union_updated[new_key] = adata

# Update the dictionary
adata_union = adata_union_updated

# Print Cells x Genes values before filtering
for name, adata in adata_union.items():
     print(f"Cells x Genes retained in {name} before filtering:  {adata.n_obs} x {adata.n_vars}")
    
print("###########################################################")
print("Filtering")
print("###########################################################")

for p in adata_union.items():
    sc.pp.filter_cells(p[1], min_genes=100)
    sc.pp.filter_genes(p[1], min_cells=3)

for name,adata in adata_union.items():
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True)

    # ribosomal genes
    pattern = r'^RP[SL]'                        
    adata.var["ribo"] = adata.var_names.str.match(pattern) 
    sc.pp.calculate_qc_metrics(adata, qc_vars=["ribo"], percent_top=[20], log1p=True, inplace=True)
                               
    #calculate the % of erithrocyte (hemoglobin)
    adata.var["Erythroid"] =  adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["Erythroid"], percent_top=[20], log1p=True, inplace=True)
    adata_union[name] = adata

np.save('./adata_union_choose_filter2.npy', adata_union) 

# MAD computation and outlier definition (genral purpose function for all continuous metrics)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric] 
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

for name, adata in adata_union.items():
    print("\n------------")
    print(name)
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 4)
        | is_outlier(adata, "log1p_n_genes_by_counts", 3)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 3)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 15#8
    )
    adata.obs["erythroid_outlier"] = adata.obs.pct_counts_Erythroid > 5
    adata.obs["n_genes_by_counts_outlier"] = adata.obs.n_genes_by_counts < 500
    
    print("low counts cells ", adata.obs.outlier.value_counts())
    print("High MT cells ", adata.obs.mt_outlier.value_counts())
    print("Erythroid cells ", adata.obs.erythroid_outlier.value_counts())
    print(f"Total number of cells for {name}: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier) & (~adata.obs.erythroid_outlier) & (~adata.obs.n_genes_by_counts_outlier) ].copy()
    if adata.n_obs > 100:
        adata_union[name] = adata
        print(f"Number of cells after filtering of low quality cells for {name}: {adata.n_obs}")
        print("------------\n")
    else:
        #Discard the adata if the object has < 100 cells
        print(f"{name} has been excluded due to low number of cells: {adata.n_obs}")
        print("------------\n")

print("###########################################################")
print("Doublet selection")
print("###########################################################")
ds_tot_py = []
dc_tot_py = []

r('''
library(Seurat)
library(scater)
library(scDblFinder)
library(BiocParallel)
''')

for name,adata in adata_union.items(): 
    data_mat = adata.X.T
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    robjects.globalenv["data_mat"] = data_mat
    r('''
        print(paste0("inizio scDblFinder: ", Sys.time()))
        print(dim(data_mat))
        sce <- scDblFinder(SingleCellExperiment(list(counts=data_mat)))
        print(paste0("fine scDblFinder: ", Sys.time())) 
        #doublet_score <- list(sce$scDblFinder.score)
        #doublet_class <- list(sce$scDblFinder.class)
        doublet_score <- sce$scDblFinder.score
        doublet_class <- sce$scDblFinder.class
    ''')
    doublet_score = list(robjects.r("doublet_score"))
    doublet_class = list(robjects.r("doublet_class"))
    ds_tot_py.append(doublet_score)
    dc_tot_py.append(doublet_class)

# Assign scDBLFinder metadata
for i,(name, adata) in enumerate(adata_union.items()):
    adata.obs["scDblFinder_score"] = ds_tot_py[i]
    adata.obs["scDblFinder_class"] = dc_tot_py[i]

for name,adata in adata_union.items():
    # Retain only singlets
    adata = adata[adata.obs['scDblFinder_class'] == "singlet"]
    adata_union[name] = adata 

print("###########################################################")
print("Combine")
print("###########################################################") 
# Combine the adata objects into one 
adata_combined = sc.concat(adata_union, join="outer")
adata_combined.write_h5ad(
  "./AdataFiltered",
  compression=hdf5plugin.FILTERS["zstd"])