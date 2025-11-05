# scRNAseq-comprehensive-analysis-pipeline
This repository aims to contain and share a series of scripts and notebooks that can be useful throughout the course of a single cell RNA sequencing analysis.
All of the following python/R scripts/notebooks are to be considered as general guidelines, sometimes the used parameters will need some fine tuning.<br />
Processes such as manual curation for Cell type identification will need human supervision.<br />
<br />
Here  below is a brief description for each script/notebook

## 0.Preprocessing,py

### First line QC
This script aimns to automate the first steps of the QC process for scRNAseq analysis.<br />
In here we extract the cont matrix starting from the filtered_feature_bc_matrix, obtained via the cellranger software (10x Genomics) and perform the first steps of quality control.<br />
To avoid the inclusion of low quality/damaged cells we remove cells according to gene content, mithocondrial gene content, erhitroyd gene content and we filter genes according to their overall expression over the entire sample.<br />
To avoid the inclusion of doublets we used "scDblFinder", a package that automatically detects doublets in the samples.
<br />
Here is briefly explained the rationale. <br />
#### Per-cell QC metrics (descriptive names):
* Library size (UMI counts per cell): total number of captured transcripts per cell.
* Molecular complexity (genes detected per cell): number of genes with at least one UMI.
* Expression concentration (top-20 gene fraction): % of total counts contributed by the 20 most expressed 
* Mitochondrial content: % of UMI mapped on mithocondrial genes (stress/apoptosis).
* Erythroid signature: % of UMI mapped on erithorid singature (contamination).
* Doublet class: output of scDblFinder identifying cells either as "doublet" or "singlet"

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
