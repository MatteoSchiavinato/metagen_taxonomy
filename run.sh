#!/usr/bin/env sh

PROJECT="M00003_Toasted_soy_diversity"
WD="/gpfs/data/fs71579/schmat90/CF/projects/meta/${PROJECT}"

echo """\
#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J TaxaDiv
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0096
#SBATCH --partition mem_0096

cd ${WD}/scripts

module load gcc

nextflow \
run \
main.nf \
-resume \
-work-dir ${WD}/scripts/work \
-with-report ${WD}/scripts/cmd.sbatch.report.html \
-with-timeline ${WD}/scripts/cmd.sbatch.timeline.html \
-with-dag ${WD}/scripts/cmd.sbatch.dag.png \
--output_dir ${WD} \
--threads 48 \
--fastq_dir ${WD}/raw_data/FINAL_trimmed_reads \
--host_genome ${WD}/raw_data/genomes/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--human_genome ${WD}/raw_data/genomes/hg38.fa \
--feed_genome ${WD}/raw_data/genomes/GCF_902167145.1.fa \
--kraken_db /binfl/lv71579/schmat90/software/kraken2/2.1.1/scripts/16S_db_silva \
--kraken_db_read_len 100 \
--min_confidence 0.25 \
--min_hit_groups 2 \
--min_counts 10 \
--max_species 10 \
--evenness shannon \
--dissimilarity bray \

""" \
> cmd.sbatch

sbatch --output cmd.sbatch.out cmd.sbatch
