# scRNAseq-comprehensive-analysis-pipeline
This repository aims to contain and share a series of scripts and notebooks that can be useful throughout the course of a single cell RNA sequencing analysis.
All of the following python/R scripts/notebooks are to be considered as general guidelines, sometimes the used parameters will need some fine tuning.
Processes such as manual curation for Cell type identification will need human supervision.__
Here  below is a brief description for each script/notebook

## 0.Preprocessing,py

This script aimns to automate the first steps of the QC process for scRNAseq analysis.__
In here we extract the cont matrix starting from the filtered_feature_bc_matrix, obtained via the cellranger software (10x Genomics) and perform the first steps of quality control.__
We filter cells according to gene content, mithocondrial gene content, erhitroyd gene content and doublet score ; we filter genes according to their overall expression over

