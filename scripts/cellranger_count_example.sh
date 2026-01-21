#!/bin/bash

#SBATCH --job-name=cellranger_count
#SBATCH --account=<your_username>
#SBATCH --mem=64GB
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --partition=ricerca
#SBATCH --error="cellranger_count.err"
#SBATCH --output="cellranger_count.out"
#SBATCH --mail-user=<your@email>
#SBATCH --mail-type=ALL

samples=(Sample1 Sample2 Sample3)

for sample in ${samples[@]}; do

mkdir ${sample}_count

/share/data/apps/cellranger/cellranger-7.1.0/cellranger count \
--id=${sample}_count \
--fastqs= <path_to_fastqs> \
--sample=${sample} \
--transcriptome=<path_to_genome_reference> \
--localcores=16

done
