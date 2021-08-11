#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J TaxaDiv
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0384
#SBATCH --partition mem_0384

cd /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/scripts

export OMP_NUM_THREADS=4
module load gcc

nextflow run main.nf -resume -work-dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/scripts/work -with-report /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/scripts/cmd.sbatch.report.html -with-timeline /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/scripts/cmd.sbatch.timeline.html -with-dag /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/scripts/cmd.sbatch.dag.png --output_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity --threads 48 --fastq_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/raw_data/FINAL_trimmed_reads --host_genome /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/raw_data/genomes/Gallus_gallus.GRCg6a.dna.toplevel.fa --human_genome /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/raw_data/genomes/hg38.fa --feed_genome /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/raw_data/genomes/GCF_902167145.1.fa --kraken_db /binfl/lv71579/schmat90/software/kraken2/2.1.1/scripts/16S_db_silva --kraken_db_read_len 100 --min_confidence 0.25 --min_hit_groups 2 --min_counts 10 --max_species 10 --evenness shannon --dissimilarity bray 

