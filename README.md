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

Another way is to build the database yourself. Maxikraken often doesn't work at the bracken stage because they don't provide the `kmer_distrib` file that is needed for bracken to work. Hence, you can build your own database and then use it. To do so, you can use `kraken2-build` in the folder where to create the database, making sure that you download both the library and the taxonomy files.

```
kraken2-build \
--download-library nt \
--download-taxonomy \
--build \
--db <your_database_path_and_name> \
--threads 16 \
--use-ftp \
```

Once this has finished, you can also run `bracken-build` in the same directory to create kmer abundance estimations.

```
bracken-build \
-k KMER_LEN \
-l READ_LEN \
-d MY_DB \
-x K_INSTALLATION \
-t THREADS
```

The `-k` option defines the k-mer length and is usually at least 1/2 of the read length. The `-l` option defines the average read length of the reads you're planning to map. The `-d` option is the path to the database that you created with `kraken2-build` (`--db` flag). The `-x` option is the path to the directory containing the `kraken2` and the `kraken2-build` executables. The `-t` option defines threads.

Once this command has finished, you have built the database properly and can use the pipeline.  


### Running the pipeline

The pipeline can be run using the `run.sh` script contained in the repository, just make sure to adjust the paths. Also, make sure you adjust the content of the `nextflow.config` file, changing the executables according to your own environment / `$PATH`. Check the default values of each command-line parameter of the pipeline and change them accordingly if you want to use something different.

##### Dependencies

This pipeline depends on the following programs in the `$PATH`. You can also not put them in the `$PATH` but in that case you have to explicitly declare their path in the `nextflow.config` file.

| Program | Version | Type          | Link                                               |
|---------|---------|---------------|----------------------------------------------------|
| Kraken2 | 2.1.1   | Program       | https://github.com/DerrickWood/kraken2/wiki/Manual |
| Bracken | 2.6.2   | Program       | https://ccb.jhu.edu/software/bracken/              |
| Python  | 3.6.*   | Interpreter   | https://www.python.org/downloads/                  |
| Pandas  | 1.0.*   | Python module | https://pandas.pydata.org/                         |
| Rscript | 3.5.*   | Interpreter   | https://cran.r-project.org/                        |


##### Run

Then, run the pipeline with a command like this one, which will also generate a report, a timeline and a direct acyclic graph for the proceeding of the pipeline itself:

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
--kraken_db /path/to/kraken_database \
--kraken_db_read_len 100 \
--min_confidence 0.25 \
--min_hit_groups 2 \
--min_counts 10 \
--max_species 10 \
```

### Output

The pipeline creates a directory called `taxonomy` inside the declared `--ouptut_dir`. Within this directory, there will be four subdirectories: `kraken`, `bracken`, `rel_abundance` and `diversity`. The first one will contain one subdirectory per analyzed sample. Each sample directory will contain the read counts per taxon as obtained by **kraken2**, which are however only intermediate files. In fact, the pipeline passes these files to **bracken** for abundance re-estimation. The results of bracken are contained in the `bracken` subdirectory.

Inside the `bracken` subdirectory, there will be one folder for each sample. Each of these sample folders will contain seven other directories: S, G, F, C, O, P, D. These stand for species, genus, family, class, order, phylum, domain. The abundance re-estimation is performed at each of these levels, so you can choose what level you prefer to consider for your analysis later on.

The files contained inside `rel_abundance` will be **relative abundance** plots (both in svg and png) based on the top 10 most represented taxa. There will be one plot for each taxonomic level (S,G,F,C,O,P,D). The associated relative abundance values are contained in the `*.rel_abundance.*` tables. The raw counts from which the relative abundance was estimated are also contained, in files that carry the keyword `*.counts.*`.

The files inside `diversity` represents diversity plots and tables [ ... to be continued ... ]
