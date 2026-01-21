# scRNAseq-comprehensive-analysis-pipeline
This project aims to collect and share a series of scripts and notebooks that can be applied throughout a single-cell RNA sequencing analysis using 10x Genomics [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) and the [scanpy](https://github.com/scverse/scanpy) python library.<br />
Note that this repository will not contain all the detailed explanation of the individual analytical processes and for this reason, it will only serve to share/display the pipeline, therefore a basic knowledge of single-cell RNA sequencing analysis is required to better understand the pipeline.
All of the following scripts/notebooks are to be considered as general guidelines, sometimes used parameters may require fine-tuning depending on the dataset.<br /> Processes such as manual curation for cell type identification will need human supervision.<br />
<br />
Below is a concise overview of each step of the workflow, with the corresponding scripts and notebooks provided.

## Data generation and Preprocessing

This section describes the upstream steps required to generate and preprocess single-cell RNA sequencing data prior to downstream analysis.  
Reference scripts and examples are provided to illustrate standard workflows based on 10x Genomics Cell Ranger, which can be integrated with the downstream analyses described in the later sections.

### FASTQ generation

Raw sequencing data in BCL format are converted to FASTQ files using standard demultiplexing workflows. In the context of 10x Genomics data, this step is typically performed using:
* [`cellranger mkfastq`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-mkfastq), a 10x Genomics wrapper around Illumina demultiplexing tools  
* Alternatively, Illumina [`bcl2fastq`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing)

This repository provides reference SLURM-based job scripts illustrating a typical HPC execution of `cellranger mkfastq`, together with a minimal example of a compatible sample sheet:
* [`cellranger_mkfastq_example.sh`](https://github.com/TomTore/scRNAseq-comprehensive-analysis-pipeline/blob/main/scripts/cellranger_mkfastq_example.sh)
* [`samplesheet_example.csv`](https://github.com/TomTore/scRNAseq-comprehensive-analysis-pipeline/blob/main/templates/samplesheet_example.csv)

Note that sample index sequences (e.g. 10x Genomics SI index sets) should always be retrieved from the official [10x Genomics documentation](https://www.10xgenomics.com/support).

### Count matrix generation and preprocessing

Gene–cell count matrices are generated from FASTQ files using the 10x Genomics alignment and quantification pipeline:
* [`cellranger count`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count#)
This step produces both `raw_feature_bc_matrix` and `filtered_feature_bc_matrix` outputs.  
The `filtered_feature_bc_matrix` is used as the primary input for all downstream analyses implemented in this repository, including quality control, integration, clustering, and annotation.
A reference SLURM-based job script illustrating a typical HPC execution of `cellranger count` is provided:
* [`cellranger_count_example.sh`](https://github.com/TomTore/scRNAseq-comprehensive-analysis-pipeline/blob/main/scripts/cellranger_count_example.sh)
All upstream scripts are provided as **reference templates**. Resource requests, module loading, file paths, and scheduler directives reflect a specific institutional HPC configuration and must be adapted to the local computing environment.
<br />

## Processing

The ["0.Processing"](https://github.com/TomTore/scRNAseq-comprehensive-analysis-pipeline/blob/main/scripts/0.Processing.py) script aims to automate the first steps of the QC process for scRNAseq analysis.<br />
In here we extract the count matrix starting from the filtered_feature_bc_matrix, obtained via the cellranger software (10x Genomics) and perform the first steps of quality control.<br />
To avoid the inclusion of low quality/damaged cells we remove cells according to gene content, mitochondrial gene content, erythroid gene content and we filter genes according to their overall expression over the entire sample.<br />
To avoid the inclusion of doublets we use [scDblFinder](https://github.com/plger/scDblFinder), a package that automatically detects doublets in the samples.
<br />
Here is briefly explained the rationale. <br />
### Per-cell QC metrics (descriptive names):
* Library size (UMI counts per cell): total number of captured transcripts per cell.
* Molecular complexity (genes detected per cell): number of genes with at least one transcript.
* Expression concentration (top-20 gene fraction): % of total counts contributed by the 20 most expressed 
* Mitochondrial content: % of UMI mapped to mitochondrial genes (stress/apoptosis).
* Erythroid signature: % of UMI mapped to erythroid genes (erythroid contamination).
* Doublet class: output of [scDblFinder](https://github.com/plger/scDblFinder) identifying cells either as "doublet" or "singlet"

### Outlier detection (robust to skew): | 𝒙_𝒊  − 𝒎𝒆𝒅(𝑿)|>𝒌∗𝑴𝑨𝑫(𝑿)

* X = Library size: k=4
* X = Molecular complexity: k=3
* X = Expression concentration (top-20): k=3
* X = Mitochondrial content: k=3 or absolute > 15%

### Hard thresholds: 𝒙_𝒊 > 𝒌

* X = Erythroid signature: k = 5%
* X = Genes detected per cell: k = 500

### scDblFinder thresholds:
* Doublet class = "doublet"

### Filtering logic:
Exclude cells flagged by any criterion above. <br />
Retain samples only if > 100 cells remain after QC.<br />
Notes: Thresholds are tunable (e.g., use higher values of k in neutrophil-rich or noisy datasets).

## From Integration to Annotation

The ["1.Integration & Annotation"](https://github.com/TomTore/scRNAseq-comprehensive-analysis-pipeline/blob/main/scripts/1.Integration%20%26%20Annotation.ipynb) notebook aims to:

1. Normalize, scale, regress, create UMAP
2. Perform data integration on multiple samples (if necessary)
3. Evaluate batch effects and QC metrics generated in the first script "0.Preprocessing.py"
4. Annotate putative cell types through a manual curation supported by an automated annotation algorithm

### Normalization

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

### Integration

Integration is the process that aims to remove unwanted batch effects arising from different samples in different biological and/or technical conditions.<br />
Data integration is advisable but it's not always required (there could be no batch effect to correct).
To perform integration we use the [harmonypy](https://github.com/slowkow/harmonypy) package, a port of the [harmony](https://github.com/immunogenomics/harmony) R package, via the custom function `do_harmony` as follows:<br />

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

### QC Metrics evaluation

To evaluate the QC metrics of the first script we now simply use the `scanpy.pl.umap` or the `scanpy.pl.embedding` (for `scanpy.pl.embedding` it is necessary to specify the `basis` parameter either as `X_umap` or as `X_harmony_umap_{batch_column}` ).<br />
If one or more of the QC metrics is all clustered together on the UMAP  (eg: mitochondrial gene content) then it would be advisable to rerun the 0.Preprocessing.py and/or the `process` and `do_harmony` functions with different parameter.

### Manual & Automated Annotation

To annotate putative cell types we use:
* Canonical Markers (i.e. cell type markers known in literature)
* A logistic regression–based classifier provided in the [celltypist](https://github.com/Teichlab/celltypist) package
* Differential gene expression analysis on the cells cluster through `scanpy.tl.rank_genes_groups`. <br />
Once all results have been interpreted, clusters can be assigned to specific cell identities, while unassigned cells may be optionally excluded or relabeled. <br />
Note that this step requires substantial time and manual effort. Be prepared to rerun clustering multiple times using different parameters, and keep in mind that the same cell type identity can be assigned to multiple clusters (thus, starting from a high-resolution clustering is often preferable).
---

Overall, the proposed framework is intended as a flexible and modular solution rather than a rigid, one-size-fits-all pipeline.  
Users are encouraged to adapt parameters, thresholds, and analytical choices to the biological context and characteristics of their own datasets. <br />

**Software versions used for testing:** Cell Ranger v7.1.0; Python 3.12.7; Scanpy v1.10.1; harmonypy v0.0.10; CellTypist v1.6.3; scDblFinder (R) v1.24.0.





