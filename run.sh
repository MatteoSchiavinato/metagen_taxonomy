#!/usr/bin/env sh

PROJECT="M00001_PacBio_chicken_gut_metagenomics"
WD="/gpfs/data/fs71579/schmat90/CF/projects/meta/${PROJECT}"

echo """\
#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J mp-TAXA
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0384
#SBATCH --partition mem_0384

cd ${WD}/scripts/wf-taxa

nextflow \
run \
main.nf \
-resume \
-work-dir ${WD}/scripts/wf-taxa/work \
-with-report ${WD}/scripts/wf-taxa/cmd.sbatch.report.html \
-with-timeline ${WD}/scripts/wf-taxa/cmd.sbatch.timeline.html \
-with-dag ${WD}/scripts/wf-taxa/cmd.sbatch.dag.png \
--output_dir ${WD} \
--threads 48 \
--ccs_dir ${WD}/filt_reads \
--extension no_host.no_human.fasta \
--kraken_db /binfl/lv71579/schmat90/software/kraken2/maxikraken2_1903_140GB \
--kraken_db_read_len 100 \
--min_confidence 0.25 \
--min_hit_groups 2 \
--min_counts 10 \
--max_species 10 \

""" \
> cmd.sbatch

sbatch --output cmd.sbatch.out cmd.sbatch
