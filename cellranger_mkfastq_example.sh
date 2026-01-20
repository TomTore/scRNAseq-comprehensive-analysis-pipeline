#!/bin/bash

#SBATCH --job-name=cellranger_mkfastq
#SBATCH --account=<your_username>
#SBATCH --mem=64GB
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --partition=ricerca
#SBATCH --error="cellranger_mkfastq.err"
#SBATCH --output="cellranger_mkfastq.out"
#SBATCH --mail-user=<your@email>
#SBATCH --mail-type=ALL

source /share/data/apps/anaconda3/bin/activate tools

/share/data/apps/cellranger/cellranger-7.1.0/cellranger mkfastq \
--id=<Folder_that_will_contain_FASTQs> \
--run=<Path_to_sequencing_run> \
--samplesheet=<samplesheet.csv> \
-p 16
