# scRNAseq-comprehensive-analysis-pipeline
This repository aims to contain and share a series of scripts and notebooks that can be useful throughout the course of a single cell RNA sequencing analysis using [scanpy](https://github.com/scverse/scanpy) python library.<br />
Note that this repository will not contain all the detailed explanation of the individual analytical processes and for this reason, it will only serve to share/display the pipeline, therefore a basic knowledge of single cell rna sequencing analysis is required to better understand the pipeline.
All of the following python/R scripts/notebooks are to be considered as general guidelines, sometimes the used parameters will need some fine tuning.<br />
Processes such as manual curation for Cell type identification will need human supervision.<br />
<br />
Here  below is a brief description for each script/notebook

## 0.Preprocessing.py

### First line QC
This script aimns to automate the first steps of the QC process for scRNAseq analysis.<br />
In here we extract the cont matrix starting from the filtered_feature_bc_matrix, obtained via the cellranger software (10x Genomics) and perform the first steps of quality control.<br />
To avoid the inclusion of low quality/damaged cells we remove cells according to gene content, mithocondrial gene content, erhitroyd gene content and we filter genes according to their overall expression over the entire sample.<br />
To avoid the inclusion of doublets we used [scDblFinder](https://github.com/plger/scDblFinder), a package that automatically detects doublets in the samples.
<br />
Here is briefly explained the rationale. <br />
#### Per-cell QC metrics (descriptive names):
* Library size (UMI counts per cell): total number of captured transcripts per cell.
* Molecular complexity (genes detected per cell): number of genes with at least one UMI.
* Expression concentration (top-20 gene fraction): % of total counts contributed by the 20 most expressed 
* Mitochondrial content: % of UMI mapped on mithocondrial genes (stress/apoptosis).
* Erythroid signature: % of UMI mapped on erithorid singature (contamination).
* Doublet class: output of [scDblFinder](https://github.com/plger/scDblFinder) identifying cells either as "doublet" or "singlet"

#### Outlier detection (robust to skew):A cell is flagged as an outlier when: ∣x−median(x)∣> k×MAD(x)
* Library size: k=4
* Molecular complexity: k=3
* Expression concentration (top-20): k=3
* Mitochondrial content: k=3 or absolute > 15%

#### Hard thresholds:
* Erythroid signature > 5%
* Genes detected per cell < 500

#### scDblFinder thresholds:
* Doublet class = "doublet"

#### Filtering logic:
Exclude cells flagged by any criterion above. <br />
Retain samples only if > 100 cells remain after QC.<br />
Notes: Thresholds are tunable (e.g., use more conservative k=5 in neutrophil-rich or noisy datasets).

## 1.Preprocessing.ipynb

### Second line QC & Annotation
This notebook aims to:

1. Normalize, scale, regress, create UMAP
2. perform data integration on multiple samples (if necessary)
3. evaluate batch effects and QC metrics generated in the first script "0.Preprocessing.py"
4. Annotate putative cell types through a manual curation supported by an automated annotation algorithm

#### Normalization

To normalize, scale, perform the regression and create the UMAP, we use a custom function called `process`. <br />
The `process` function calls several other modules and functions provided by the scanpy package as follows: <br />

```
def process(adata, resolution=0.4):    
    # Check if data is raw
    if adata.X.toarray().max() > 100:
        #save raw data to "raw_counts"
        adata.X = adata.layers["raw_counts"].copy()
        # Normalize each cell by total counts over all genes
        sc.pp.normalize_total(adata, target_sum=10**4)
        # Logarithmize the data matrix
        sc.pp.log1p(adata)
    else:
        print("Data are already log trasformed")
    # Annotate highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if adata.raw is None:
        print("log1p data before filtering have been saved in .raw.X")
        adata.raw = adata
    # Save norm data into log1p
    adata.layers["log1p"] = adata.X.copy()
    # Restrict to HVGs
    adata = adata[:, adata.var['highly_variable']]
    # Uses simple linear regression to "regress out" (mostly) unwanted sources of variation
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # Scale data to unit variance and zero mean
    sc.pp.scale(adata, max_value=10)
    #  Computes Principal Component Analysis (PCA) coordinates, loadings and variance decomposition.
    sc.tl.pca(adata, svd_solver='arpack', random_state=1)
    # Compute the nearest neighbors distance matrix and a neighborhood graph of observations.
    sc.pp.neighbors(adata, random_state=1)
    # Cluster cells using the Leiden algorithm
    sc.tl.leiden(adata, resolution=resolution, random_state=1, key_added= f"leiden_{resolution}")
    # Embed the neighborhood graph using Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction
    sc.tl.umap(adata, random_state=1)
    
    return adata
```
<br />
#### Integration
Integration is the process that aims to remove unwanted batch effects arising from different samples in different biological and/or technical conditions.<br />
Data integration is advisable but it's not always necessary (there could be no batch effect to correct).
To perform integration we use the [harmonypy package](https://github.com/slowkow/harmonypy) (a port of the [harmony](https://github.com/immunogenomics/harmony) R package ) via the custom function `do_harmony` as follows:<br />

```

def do_harmony(adata_to_harmonize, resolution=0.4, batch_column='orig.ident'):
    adata = adata_to_harmonize.copy()
    # Create "batch" as the batch contianing column to use for harmony
    adata.obs['batch'] = adata.obs[batch_column].astype(str)
    # Save original UMAP
    adata.uns['orig_umap'] = adata.uns['umap']
    adata.obsm['orig_X_umap'] = adata.obsm['X_umap']

    # running harmony
    ho = hm.run_harmony(
        adata.obsm['X_pca'],     # embedding PCA (n_cells x n_pcs)
        adata.obs,               # DataFrame con colonna batch
        'batch', random_state=42  # nome della colonna batch
    )
    # Traspose to have shape (n_cells x n_pcs)
    adata.obsm[f'X_harmony_{batch_column}'] = ho.Z_corr.T

    # Generate neighbors and UMAP with the harmony correction
    sc.pp.neighbors(adata, random_state=42, use_rep=f'X_harmony_{batch_column}', key_added=f'harmony_nn_{batch_column}')
    sc.tl.leiden(adata, random_state=42, resolution=resolution, neighbors_key=f'harmony_nn_{batch_column}', key_added=f"harmony_leiden_{batch_column}_{resolution}")
    sc.tl.umap(adata, random_state=42, neighbors_key=f'harmony_nn_{batch_column}')

    # Save UMAP on X_harmony_umap
    adata.uns[f'harmony_umap_{batch_column}'] = adata.uns['umap']
    adata.obsm[f'X_harmony_umap_{batch_column}'] = adata.obsm['X_umap']

    # Set the first umap (not batch corrected) on the default location
    adata.uns['umap'] = adata.uns['orig_umap']
    adata.obsm['X_umap'] = adata.obsm['orig_X_umap']
    del adata.uns['orig_umap']
    del adata.obsm['orig_X_umap']

    return adata
```

#### QC Metrics evaluation

To evaluate the QC metrics of the first script we now simply use the `scanpy.pl.umap` or the `scanpy.pl.embedding` (for `scanpy.pl.embedding` it is necessary to specify the `basis` parameter either as `X_umap` or as `X_harmony_umap_{batch_column}` ).<br />
If one or more of the QC metrics is all clustered together on the UMAP  (eg: mithocondrial gene content) then it would be advisable to rerun the 0.Preprocessing.py and/or the `process` and `do_harmony` functions with different parameter.

#### Manual & Automated Annotation

To annotate putative cell types we use:
* canonical Markers (i.e. cell type markers known in literature)
* an automated annotation logistic classifier provided in the [celltypist](https://github.com/Teichlab/celltypist) package
* Differential gene expression analysis on the cells cluster through `scanpy.tl.rank_genes_groups`. <br />
Once all the results have been interpreted cell identities can be correctly assigned to the clusters and eliminate the unknown cells. <br />
Note that this step requires a lot of time and dedication, be prepared to rerun clustering multiple times with different parameters and remember that you can assign the same cell type identity to multiple clusters (so it's usually better to start with an high resolution clustering)





