#!/usr/bin/env sh

WD="/gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity"

# ------------------------------------------------------------------------------
# count raw reads
cd $WD/raw_data/FINAL_trimmed_reads

{ \
for i in `ls | egrep "fastq$"`
do
  COUNT=$(bioawk -c fastx 'END {print NR}' $i); NAME=$(echo $i | sed -e 's/.fastq//')
  echo -e "$NAME\t$COUNT"
done; } | \
sort \
> read_counts.tsv

# ------------------------------------------------------------------------------
# count filtered reads
cd $WD/filt_reads

{ \
for i in `ls | egrep "fasta$"`
do
  COUNT=$(bioawk -c fastx 'END {print NR}' $i); NAME=$(echo $i | sed -e 's/.fasta//')
  echo -e "$NAME\t$COUNT"
done; } | \
sort \
> read_counts.tsv

# ------------------------------------------------------------------------------
# combine sample information and sample diversity values
cd $WD/diversity
if [ ! -d plots ]; then mkdir plots; fi
cd plots

python3.8 \
${WD}/scripts/src/combine-information.py

# ------------------------------------------------------------------------------
# plot diversity

# cd $WD/diversity/plots

cd /gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/TEST/plots

Rscript \
${WD}/scripts/src/plot-diversity.Rscript \
RES.G.counts.diversity.tsv \
. \


# alpha diversity

### RICHNESS + EVENNESS || CECUM || VARIETY ###
### RICHNESS + EVENNESS || ILEUM || VARIETY ###
### RICHNESS + EVENNESS || CECUM || TEMPERATURE ###
### RICHNESS + EVENNESS || ILEUM || TEMPERATURE ###
### RICHNESS + EVENNESS || CECUM vs ILEUM ###

# beta diversity

### DISSIMILARITY || CECUM || VARIETY ###
### DISSIMILARITY || ILEUM || VARIETY ###
### DISSIMILARITY || CECUM || TEMPERATURE ###
### DISSIMILARITY || ILEUM || TEMPERATURE ###
### DISSIMILARITY || CECUM vs ILEUM ###

# statistical testing

### KRUSKAL-WALLIS + MANN WHITNEY on SEX in CECUM || VARIETY ###
### KRUSKAL-WALLIS + MANN WHITNEY on SEX in CECUM || TEMPERATURE ###
### KRUSKAL-WALLIS + MANN WHITNEY on SEX in ILEUM || VARIETY ###
### KRUSKAL-WALLIS + MANN WHITNEY on SEX in ILEUM || TEMPERATURE ###
### KRUSKAL-WALLIS + MANN WHITNEY on CECUM vs ILEUM without SEX ###
