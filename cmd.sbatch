#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J mp-TAXA
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0384
#SBATCH --partition mem_0384

cd /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-taxa

nextflow run main.nf -resume -work-dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-taxa/work -with-report /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-taxa/cmd.sbatch.report.html -with-timeline /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-taxa/cmd.sbatch.timeline.html -with-dag /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-taxa/cmd.sbatch.dag.png --output_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics --threads 48 --ccs_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/filt_reads --extension no_host.no_human.fasta --kraken_db /gpfs/data/fs71579/schmat90/chicken/kraken2/KRAKEN_DB --kraken_db_read_len 100 --min_confidence 0.25 --min_hit_groups 2 --min_counts 10 --max_species 10 

