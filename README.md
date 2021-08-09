# Pipeline for the taxonomy analysis of a set of metagenomic reads

This pipeline takes as input a set of one or more read files (single-end illumina, pacbio) and performs a series of analysis with them, aimed at characterizing the taxa contained in the read set. The taxonomic classification tools used are *alignment free*, hence to use paired end reads just pass them as independent files.

Short reads are usually quality-trimmed removing adapter content and low-quality reads. The process (trimming) is usually done with tools like [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) or [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Make sure you trim your reads with these tools before you use this pipeline.

Metagenomic reads are usually cleaned removing any read mapping against the host genome (e.g. *Gallus gallus*), against the human genome (*Homo sapiens*), and against the main source of food of the host (e.g. *Zea mays*). The pipeline first maps the reads against a concatenated version of these genomes, and the unmapped reads are extracted with `samtools view -b -h -f 0x4`. The `-f 0x4` flag indicates to samtools that it has to retain only reads that did **not** map against the reference. The resulting **bam** file is then passed to `samtools fasta` to obtain a fasta file. The fasta files obtained from this filtering step are those that are then passed to the taxonomic classification part of the pipeline.

The pipeline assigns the reads to taxa using [Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual), and then refines read counts statistically using [Bracken](https://github.com/jenniferlu717/BrackenA). The refined counts are used to estimate richness, diversity and relative abundance using custom **python3** and **Rscript** scripts, which can be found inside the `/src` directory and are executed by the pipeline automatically.

### Preparing the database for Kraken2

Kraken2 depends on a database that has to be pre-built. Depending on your needs, there are many databases that you can use. Here are the three main databases.

##### Traditional database

Kraken2 requires a database in order to classify the reads. This pipeline assumes you have it.

The classic usage requires you to build the database yourself. To do so, you can use `kraken2-build` in the folder where to create the database, making sure that you download both the library and the taxonomy files.

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

##### 16S databases

On of the most common usages of the kraken2+bracken pipeline is to detect enrichment of 16S ribosomal RNAs. For this purpose, one should not use a traditional database built from NCBI, but rather any of the three main 16S rRNA databases: [Silva](https://www.arb-silva.de/), [Greengenes](https://greengenes.secondgenome.com/), or [RDP](https://rdp.cme.msu.edu/).

These three databases are provided inside the installation of kraken2. In the **kraken2** main directory there is a subdirectory called `scripts` which contains many useful scripts to build and prepare databases (and other things). Depending on which of the three databases you want to use, you can run any of the corresponding script (which are called, for example, `16S_silva_installation.sh`).

These scripts, however, are not perfectly fit for running them from the shell "as is". In fact, they use two variables called `$KRAKEN2_DB_NAME` and `$KRAKEN2_THREAD_CT` which are not initialized in the script. These two probably come from the command line, from another calling script, but can be declared in the script too. Here's what I added to the script:

```
KRAKEN2_DB_NAME=$@
KRAKEN2_THREAD_CT=16
```

This way, all you have to do is to pass a string as *first and only* argument of the script when you run it and that will be the directory containing your database at the end. The other variable simply sets the number of threads for the building process, so I set it to 16 but it can be any number that fits your hardware.

This was the final command for the preparation of the **Silva** database:

```
sh 16S_silva_installation.sh 16S_db_silva
cd 16S_db_silva
bracken-build -k 35 -l 100 -t 16 -d .
```

This command runs the building script from the kraken2 suite. Then it goes into the created directory and it builds a bracken database on top of it. The `*.sh` script provided with kraken2 builts a kraken2 database with a kmer length of 35: hence, the bracken database has to be built with the same kmer length (`-k 35`). The read length depends on your dataset. A normal illumina dataset likely contains reads that are ~ 100 nt long, so `-l 100` was set. If your reads are of a different length, make sure to change this parameter. The threads are set with `-t 16`, and the database location is set with `-d`. Finally, `bracken-build` requires the kraken2 executables in the `$PATH`. If these are **not** in the path, then you can declare a parameter called `-x` which specifies the actual full path where these executables can be found (i.e. the base directory of the kraken2 program in your filesystem).

If you prefer to use **GreenGenes** or **RDP**, simply do the same operations as I did here for **Silva** but on the other two scripts.

##### MaxiKraken

One of the many available kraken2 databases is **maxikraken2**. This database contains ~140 GB of data from metagenomic species. Be aware that Maxikraken often doesn't work at the bracken stage because they don't provide the `kmer_distrib` file that is needed for bracken to work. Hence, use maxikraken2 only if you're *not* planning to run bracken, or if you are willing to spend some time to build the bracken database yourself (it may take ages and not work swiftly, it's a buggy process).

The following commands must be executed in the directory where you want to create the database, in order to generate it:

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

The pipeline can be run using the `run.sh` script contained in the repository, just make sure to adjust the paths. Also, make sure you adjust the content of the `nextflow.config` file, changing the executables according to your own environment / `$PATH`. Check the default values of each command-line parameter of the pipeline and change them accordingly if you want to use something different.

##### Dependencies

This pipeline depends on the following programs in the `$PATH`. You can also not put them in the `$PATH` but in that case you have to explicitly declare their path in the `nextflow.config` file.

| Program  | Version | Type          | Link                                                               |
|----------|---------|---------------|--------------------------------------------------------------------|
| Nextflow | 21.04.1 | Interpreter   | https://www.nextflow.io/                                           |
| Seqtk    | 1.3     | Program       | https://anaconda.org/bioconda/seqtk                                |
| Minimap2 | 2.18    | Program       | https://github.com/lh3/minimap2                                    |
| Hisat2   | 2.1.1   | Program       | http://daehwankimlab.github.io/hisat2/download/#version-hisat2-221 |
| Samtools | 1.12    | Program       | http://www.htslib.org/download/                                    |  
| Kraken2  | 2.1.1   | Program       | https://github.com/DerrickWood/kraken2/wiki/Manual                 |
| Bracken  | 2.6.2   | Program       | https://ccb.jhu.edu/software/bracken/                              |
| Python   | 3.6.*   | Interpreter   | https://www.python.org/downloads/                                  |
| Pandas   | 1.0.*   | Python module | https://pandas.pydata.org/                                         |
| Rscript  | 3.5.*   | Interpreter   | https://cran.r-project.org/                                        |
| dplyr    | 1.0.6   | R package     | https://dplyr.tidyverse.org/                                       |
| ggplot2  | 3.3.3   | R package     | https://ggplot2.tidyverse.org/reference/                           |
| ggpubr   | 0.4.0   | R package     | https://rpkgs.datanovia.com/ggpubr/                                |
| reshape2 | 1.4.4   | R package     | https://www.rdocumentation.org/packages/reshape2/versions/1.4.4    |
| rstatix  | 0.7.0   | R package     | https://cran.r-project.org/web/packages/rstatix/index.html         |
| vegan    | 2.5.*   | R package     | https://cran.r-project.org/web/packages/vegan/vegan.pdf            |


##### Run

To see the pipeline help section, run:

```
nextflow run main.nf --help
```

Then, run the pipeline with a command like this one, which will also generate a report, a timeline and a direct acyclic graph for the proceeding of the pipeline itself:

```
nextflow \
run \
main.nf \
--output_dir /path/to/output \
--threads 48 \
--ccs_dir /path/to/reads \
--extension fasta \
--kraken_db /path/to/kraken_database \
--kraken_db_read_len 100 \
--min_confidence 0.25 \
```

##### Nextflow tips

- You can resume a crashed run by re-running the same command and adding `-resume` as an option after `run`. More info [here](https://www.nextflow.io/docs/latest/getstarted.html).
- You can specify where to save the temporary files of the pipeline by specifying a `-work-dir` directory right after `run`. More info on that and on other options available in `nextflow run` can be found [here](https://www.nextflow.io/docs/latest/cli.html#clean).
- The `work` directory tends to become quite crowded and full thousands of internal nextflow files, so every now and then make sure you clean it with `nextflow clean`. More info [here](https://www.nextflow.io/docs/latest/cli.html#clean)

##### Editing the pipeline

If you don't like the plots, if you want to change something in the code to accustom it to your own taste, you can edit directly the `Rscript` or `py` files inside the `/src` directory of the git repository that you cloned.

### Output

The pipeline creates a directory called `taxonomy` inside the declared `--ouptut_dir`. Within this directory, there will be four subdirectories: `kraken`, `bracken`, `rel_abundance` and `diversity`.

##### Read classification

The `kraken` directory will contain one subdirectory per analyzed sample. Each sample directory will contain the read counts per taxon as obtained by **kraken2**, which are however only intermediate files. In fact, the pipeline passes these files to **bracken** for abundance re-estimation. The results of bracken are contained in the `bracken` subdirectory.

Inside the `bracken` subdirectory, there will be one folder for each sample. Each of these sample folders will contain seven other directories: S, G, F, C, O, P, D. These stand for species, genus, family, class, order, phylum, domain. The abundance re-estimation is performed at each of these levels, so you can choose what level you prefer to consider for your analysis later on.

Using the `--genus_only` option, the workflow will perform a kraken2+bracken analysis only at the genus level which is the most common option. Otherwise, it will perform an analysis at any taxonomic rank (S,G,F,C,O,P,D).

##### Relative abundance

The files contained inside `rel_abundance` will be **relative abundance** tables (`*.rel_abundance.*`) containing the relative abundance values. The raw counts from which the relative abundance was estimated are also contained, in files that carry the keyword `*.counts.*`.

##### Alpha and beta diversity

The files inside `diversity` represents **diversity** tables. Each table represents a diversity analysis performed at a different taxonomic level (S,G,F,C,O,P,D). The diversity scores are obtined from the `*.counts.*` tables found inside the `rel_abundance` directory. The diversity analysis performed includes alpha diversity (richness, evenness) and beta diversity (dissimilarity). The evenness values depend on the setting chosen with `--evenness`. The setting influences the type of calculation performed, see R documentation for details. The same applies to the `--dissimilarity` option. Default settings are `--evenness shannon` and `--dissimilarity bray`, which are commonly employed in metagenomic workflows.

Using the `--skip_diversity`  option, the tool will only perform a kraken2+bracken analysis without performing any diversity analysis. This could reduce computational times but returns incomplete results. 
