# Pipeline for the taxonomy analysis of a set of metagenomic reads

This pipeline takes as input a set of one or more read files (single-end illumina, pacbio) and performs a series of analysis with them, aimed at characterizing the taxa contained in the read set. The illumina reads can for now only be **single-end** due to how the code is implemented, but it will be adjusted in the future. For now, treat the two read files as independent. The pipeline assigns the reads to taxa using [Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual), and then refines read counts statistically using [Bracken](https://github.com/jenniferlu717/BrackenA). The refined counts are used to estimate richness, diversity and relative abundance using custom **python3** and **Rscript** scripts, which can be found inside the `/src` directory and are executed by the pipeline automatically.

### Cleaning the reads

The pipeline assumes that the reads you use do not need any taxonomic cleaning or quality trimming. If you need to do so, take care to do it yourself.

Metagenomic reads are usually cleaned removing any read mapping against the host genome (e.g. *Gallus gallus*) and against the human genome. Reads are mapped against a concatenated version of these two genomes, and the unmapped reads are extracted with `samtools view -b -h -f 0x4`. The `-f 0x4` flag indicates to samtools that it has to retain only reads that did **not** map against the reference. The resulting **bam** file can then be passed to `samtools fasta` to obtain a fasta file (or any other tool that converts from bam to fasta, there are many).

Short reads are usually quality-trimmed removing adapter content and low-quality reads. The process (trimming) is usually done with tools like [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) or [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Make sure you trim your reads with these tools.

### Preparing the database for Kraken2

Kraken2 requires a database in order to classify the reads. This pipeline assumes you have it. One of the many available kraken2 databases is **maxikraken2**. This database contains ~140 GB of data from metagenomic species. The following commands must be executed in the directory where you want to create the database, in order to generate it:

```
mkdir maxikraken2_1903_140GB
cd maxikraken2_1903_140GB
wget -c https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/hash.k2d
wget https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/opts.k2d
wget https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/taxo.k2d

kraken2-inspect --db .
```

The first commands create a folder and download the three components of the database inside. The inspect command from kraken2 instead checks if everything is in place. It can take up to 30 minutes, so run it when you can have a break. It outputs a list of species found and other details about the database which may be useful to your publications, so save its output to a file.

### Running the pipeline

The pipeline can be run using the `run.sh` script contained in the repository, just make sure to adjust the paths. Also, make sure you adjust the content of the `nextflow.config` file, changing the executables according to your own environment / `$PATH`. Check the default values of each command-line parameter of the pipeline and change them accordingly if you want to use something different. Then, run the pipeline with:

```
nextflow run main.nf
```

If you don't want to change the content of the `nextflow.config` file, just pass its config parameters as command line parameters, like in this example:

```
nextflow \
run \
main.nf \
-resume \
-work-dir /path/to/work \
-with-report /path/to/cmd.sbatch.report.html \
-with-timeline /path/to/cmd.sbatch.timeline.html \
-with-dag /path/to/cmd.sbatch.dag.png \
--output_dir /path/to/output \
--threads 48 \
--ccs_dir /path/to/reads \
--extension fasta \
--kraken_db /path/to/maxikraken2_1903_140GB \
--kraken_db_read_len 100 \
--min_confidence 0.25 \
--min_hit_groups 2 \
--min_counts 10 \
--max_species 10 \
```
