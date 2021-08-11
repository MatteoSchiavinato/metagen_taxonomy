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
# create file with read counts and fractions 
# ready for publication 

cd ${WD}

{ \
echo -e "Sample\tRaw\tFiltered\tFrac"; \
paste <(cat raw_data/FINAL_trimmed_reads/read_counts.tsv | sort) \
<(cat filt_reads/read_counts.tsv | sed -e 's/.no-hhf//' | sort) | \
awk '{print $1"\t"$3"\t"$2"\t"$4}' | \
cut -f 1,3,4 | \
awk '{print $0"\t"$3/$2*100}'; } \
> manuscript/File_READS.tsv

# ------------------------------------------------------------------------------
# count annotated reads at any taxonomic rank 
cd ${WD}/taxonomy/rel_abundance

unset X
declare -a X=(D P O C F G S)
for TAXRANK in ${X[@]}
do
paste \
<(cat RES.${TAXRANK}.rel_abundance.tsv | grep "Taxa" | cut -f 2- | tr "\t" "\n") \
<(cat RES.${TAXRANK}.rel_abundance.tsv | grep "Unassigned" | cut -f 2- | tr "\t" "\n" | awk '{print 100-$1}') \
> RES.${TAXRANK}.annotated.tsv
done

cd ${WD}/taxonomy/rel_abundance

{ \
echo -e "Sample\tDomain\tPhylum\tOrder\tClass\tFamily\tGenus\tSpecies"; \
paste \
RES.D.annotated.tsv RES.P.annotated.tsv RES.O.annotated.tsv RES.C.annotated.tsv \
RES.F.annotated.tsv RES.G.annotated.tsv RES.S.annotated.tsv | \
cut -f 1,2,4,6,8,10,12,14; } \
> ../../manuscript/File_ANNOTATED.tsv

# ------------------------------------------------------------------------------
# combine sample information and sample diversity values
cd $WD/taxonomy/diversity
if [ ! -d plots ]; then mkdir plots; fi
cd plots

python3.8 \
${WD}/scripts/src/combine-information.py

# ------------------------------------------------------------------------------
# plot diversity

cd $WD/taxonomy/diversity/plots

Rscript \
${WD}/scripts/src/plot-diversity.Rscript \
RES.G.counts.diversity.tsv \
. 

cp taxonomy/diversity/RES.F.counts.diversity.tsv manuscript/File_DIVERSITY_FAMILY.tsv
cp taxonomy/diversity/RES.G.counts.diversity.tsv manuscript/File_DIVERSITY_GENUS.tsv
cp taxonomy/diversity/RES.P.counts.diversity.tsv manuscript/File_DIVERSITY_PHYLUM.tsv
