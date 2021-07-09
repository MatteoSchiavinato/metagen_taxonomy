#!/usr/bin/env nextflow

// help section
if (params.help) {

  println """

  Matteo schiavinato
  Metagenomic ccs taxa assignment
  June 2021-July 2021

  ### misc ###

  --threads               Number of parallel threads

  ### input ###

  --ccs_dir               Directory containing FASTA files with pacbio ccs reads
  --fastq_dir             Directory containing FASTQ files with short reads
  --fasta_dir             Directory containing FASTA files with short reads

  ### removal of unwanted reads ###

  --host_genome           E.g. Gallus gallus (none)
  --human_genome          E.g. Homo sapiens (none)
  --feed_genome           E.g. Zea mays (none)

  ### diversity ###

  --evenness              One of: "shannon", "simpson", "invsimpson"
                          (shannon)
  --dissimilarity         One of: "bray", "jaccard", "manhattan", "euclidean",
                          "canberra", "clark", "kulczynski", "gower", "altGower",
                          "morisita", "horn", "mountford", "raup", "binomial",
                          "chao", "cao", "mahalanobis", "chisq", "chord"
                          (bray)

  ### plotting ###

  --max_species           Max. number of species to include in the plot (10)

  ### output ###

  --output_dir            Base directory where all the sub-directories with the outputs
                          will be placed.
                          Ideally, this is the directory containing /raw_data and /scripts


  """
  exit 0
}

// -----------------------------------------------------------------------------

// read fasta files from long pacbio reads
if (params.ccs_dir) {

  Channel
  .fromPath("${params.ccs_dir}/*.{fa,fasta}")
  .map{ it -> [ it.simpleName, it] }
  .set{ Ccs }

  // remove unwanted reads
  // human

  process remove_human_ccs {

    executor = "local"
    cpus = params.threads

    input:
      tuple val(sample_id), file(fasta) from Ccs

    output:
      tuple val(sample_id), file("${sample_id}.no-h.fasta") into Ccs_in_nohuman

    script:
      """
      ${MINIMAP2} \
      -H \
      -a \
      -t ${params.threads} \
      ${params.human_genome} \
      ${fasta} | \
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      - \
      > ${sample_id}.no-h.fasta \

      """
  }

  // host

  process remove_host_ccs {

    executor = "local"
    cpus = params.threads

    input:
      tuple val(sample_id), file(fasta) from Ccs_in_nohuman

    output:
      tuple val(sample_id), file("${sample_id}.no-hh.fasta") into Ccs_in_nohuman_nohost

    script:
      """
      ${MINIMAP2} \
      -H \
      -a \
      -t ${params.threads} \
      ${params.host_genome} \
      ${fasta} | \
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      - \
      > ${sample_id}.no-hh.fasta \

      """
  }

  // feed

  process remove_feed_ccs {

    executor = "local"
    cpus = params.threads

    publishDir "${params.output_dir}/filt_reads", mode: "copy"

    input:
      tuple val(sample_id), file(fasta) from Ccs_in_nohuman_nohost

    output:
      tuple val(sample_id), file("${sample_id}.no-hhf.fasta") into Ccs_in_nohuman_nohost_nofeed

    script:
      """
      ${MINIMAP2} \
      -H \
      -a \
      -t ${params.threads} \
      ${params.feed_genome} \
      ${fasta} | \
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      - \
      > ${sample_id}.no-hhf.fasta \

      """
  }
} else {
  Ccs_in_nohuman_nohost_nofeed = Channel.empty()
}


// build hisat2 index if needed
if (params.fastq_dir || params.fasta_dir) {

  process build_hisat2_index {

    executor = "local"
    cpus = params.threads

    output:
      tuple val("human_genome"), path("human_genome*") into Human_genome_idx
      tuple val("host_genome"), path("host_genome*") into Host_genome_idx
      tuple val("feed_genome"), path("feed_genome*") into Feed_genome_idx

    script:
      """
      ${HISAT2_BUILD} -p ${params.threads} ${params.human_genome} human_genome
      ${HISAT2_BUILD} -p ${params.threads} ${params.host_genome} host_genome
      ${HISAT2_BUILD} -p ${params.threads} ${params.feed_genome} feed_genome
      """
    }
}

// read fastq files from short reads
if (params.fastq_dir) {

  Channel
    .fromPath("${params.fastq_dir}/*.{fq,fastq}")
    .map{ it -> [ it.simpleName, it] }
    .set{ Fastq }


  // convert reads to fasta format

  process convert_to_fasta {

    executor = "local"
    cpus = 1

    input:
      tuple val(sample_id), file(fastq) from Fastq

    output:
      tuple val(sample_id), file("${sample_id}.fa") into Fastq_to_fasta

    script:
      """
      ${SEQTK} \
      seq \
      -A \
      ${fastq} \
      > ${sample_id}.fa \

      """
  }
} else {
  Fastq_to_fasta = Channel.empty()
}


// read fasta files from short reads
if (params.fasta_dir) {

  Channel
  .fromPath("${params.fasta_dir}/*.{fa,fasta}")
  .map{ it -> [ it.simpleName, it] }
  .set{ Fasta }

} else {
  Fasta = Channel.empty()
}


if (params.fastq_dir || params.fasta_dir) {

  // concatenate channels

  Fastq_to_fasta
  .concat(Fasta)
  .set{ Fastq_and_fasta }

  // remove unwanted reads
  // human

  process remove_human_reads {

    executor = "local"
    cpus = params.threads

    input:
      tuple val(sample_id), file(fasta) from Fastq_and_fasta
      tuple val(index_prefix), path(index) from Human_genome_idx

    output:
      tuple val(sample_id), file("${sample_id}.no-h.fasta") into Fastq_and_fasta_in_nohuman

    script:
      """
      ${HISAT2} \
      -f \
      -x ${index_prefix} \
      -U ${fasta} \
      -S ${sample_id}.no-h.sam &&
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      ${sample_id}.no-h.sam \
      > ${sample_id}.no-h.fasta \

      """
  }

  // host

  process remove_host_reads {

    executor = "local"
    cpus = params.threads

    input:
      tuple val(sample_id), file(fasta) from Fastq_and_fasta_in_nohuman
      tuple val(index_prefix), path(index) from Host_genome_idx

    output:
      tuple val(sample_id), file("${sample_id}.no-hh.fasta") into Fastq_and_fasta_in_nohuman_nohost

    script:
      """
      ${HISAT2} \
      -f \
      -x ${index_prefix} \
      -U ${fasta} \
      -S ${sample_id}.no-hh.sam &&
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      ${sample_id}.no-hh.sam \
      > ${sample_id}.no-hh.fasta \

      """
  }

  // feed

  process remove_feed_reads {

    executor = "local"
    cpus = params.threads

    publishDir "${params.output_dir}/filt_reads", mode: "copy"

    input:
      tuple val(sample_id), file(fasta) from Fastq_and_fasta_in_nohuman_nohost
      tuple val(index_prefix), path(index) from Feed_genome_idx

    output:
      tuple val(sample_id), file("${sample_id}.no-hhf.fasta") into Fastq_and_fasta_in_nohuman_nohost_nofeed

    script:
      """
      ${HISAT2} \
      -f \
      -x ${index_prefix} \
      -U ${fasta} \
      -S ${sample_id}.no-hhf.sam &&
      ${SAMTOOLS} fasta \
      -@ ${params.threads} \
      -f 0x4 -F 0x0100 \
      ${sample_id}.no-hhf.sam \
      > ${sample_id}.no-hhf.fasta \

      """
  }
} else {
  Fastq_and_fasta_in_nohuman_nohost_nofeed = Channel.empty()
}


// mix pacbio and illumina streams
Ccs_in_nohuman_nohost_nofeed
  .concat(Fastq_and_fasta_in_nohuman_nohost_nofeed)
  .set{ Reads_in }



// kraken2

process kraken2 {

  cpus = params.threads
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/kraken2", mode: "copy"

  input:
    tuple val(sample_id), file(fasta) from Reads_in

  output:
    path "${sample_id}" into Kraken2_out
    tuple val(sample_id), file("${sample_id}/${sample_id}.report") into Counts
    tuple val(sample_id), file("${sample_id}/${sample_id}.kraken.log") into Kraken_logs

  script:
    """
    if [ ! -d ${sample_id} ]; then mkdir ${sample_id}; fi &&
    { \
    ${KRAKEN2} \
    --db ${params.kraken_db} \
    --threads ${params.threads} \
    --unclassified-out ${sample_id}/${sample_id}.uncl.fa \
    --output ${sample_id}/${sample_id}.out \
    --confidence ${params.min_confidence} \
    --report ${sample_id}/${sample_id}.report \
    --minimum-hit-groups ${params.min_hit_groups} \
    ${fasta}; } \
    &> ${sample_id}/${sample_id}.kraken.log \

    """
}

// count classified

process count_classified {

  cpus = 1
  executor = "local"

  input:
    tuple val(sample_id), file(kraken_log) from Kraken_logs

  output:
    file "${sample_id}.classif_reads.tsv" into Classif_reads

  script:
    """
    CL=\$(cat ${kraken_log} | tr -s " " | grep "sequences classified" | cut -d " " -f 2) &&
    CL_FRAC=\$(cat ${kraken_log} | tr -s " " | grep "sequences classified" | cut -d " " -f 5 | tr -d "()%") &&
    UNCL=\$(cat ${kraken_log} | tr -s " " | grep "sequences unclassified" | cut -d " " -f 2) &&
    UNCL_FRAC=\$(cat ${kraken_log} | tr -s " " | grep "sequences unclassified" | cut -d " " -f 5 | tr -d "()%") &&
    echo -e "${sample_id}\t\${CL}\t\${CL_FRAC}\t\${UNCL}\t\${UNCL_FRAC}" \
    > ${sample_id}.classif_reads.tsv

    """
}

// create classified file

process create_classified_table {

  cpus = 1
  executor = "local"
  publishDir "${params.output_dir}/taxonomy", mode: "copy"

  input:
    file classif_tables from Classif_reads.collect()

  output:
    file "classified_reads.tsv" into classified_reads

  script:
    """
    { echo -e 'Sample\tC\tC_frac\tU\tU_frac' &&
    cat ${classif_tables}; } > classified_reads.tsv
    """
}


// bracken

process bracken {

  cpus = params.threads
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/bracken", mode: "copy"

  input:
    tuple val(sample_id), file(report) from Counts

  output:
    path "${sample_id}" into Bracken_out
    tuple val("S"), file("${sample_id}/S/${sample_id}.S.bracken.report") into Bracken_S
    tuple val("G"), file("${sample_id}/G/${sample_id}.G.bracken.report") into Bracken_G
    tuple val("F"), file("${sample_id}/F/${sample_id}.F.bracken.report") into Bracken_F
    tuple val("C"), file("${sample_id}/C/${sample_id}.C.bracken.report") into Bracken_C
    tuple val("O"), file("${sample_id}/O/${sample_id}.O.bracken.report") into Bracken_O
    tuple val("P"), file("${sample_id}/P/${sample_id}.P.bracken.report") into Bracken_P
    tuple val("D"), file("${sample_id}/D/${sample_id}.D.bracken.report") into Bracken_D

  script:
    """
    if [ ! -d ${sample_id} ]; then mkdir ${sample_id}; fi &&
    unset X &&
    declare -a X=(S G F O C P D) &&
    for LEVEL in \${X[@]}
    do
      if [ ! -d ${sample_id}/\${LEVEL} ]; then mkdir ${sample_id}/\${LEVEL}; fi &&
      { ${BRACKEN} \
      -d ${params.kraken_db} \
      -i ${report} \
      -o ${sample_id}/\${LEVEL}/${sample_id}.\${LEVEL}.bracken \
      -w ${sample_id}/\${LEVEL}/${sample_id}.\${LEVEL}.bracken.report \
      -r ${params.kraken_db_read_len} \
      -l \${LEVEL} \
      -t ${params.min_counts} || \
      echo "This table wasn't generated because no reads were classified at this taxonomic level."; } \
      &> ${sample_id}/\${LEVEL}/${sample_id}.bracken.log
    done \

    """
}


// combine channels

Bracken_S
  .mix(Bracken_G, Bracken_F, Bracken_C, Bracken_O, Bracken_P, Bracken_D)
  .groupTuple()
  .map{ it -> [ it[0], it[1].collect() ] }
  .set{ Bracken_reports }


// compute relative abundance

process plot_relative_abundance {

  cpus = 1
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/rel_abundance", mode: "copy"

  input:
    tuple val(level), file(reports) from Bracken_reports
    file classified_reads

  output:
    file "RES.${level}.rel_abundance.tsv" into bracken_frac
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv" into bracken_top
    tuple val(level), file("RES.${level}.counts.tsv") into Bracken_counts
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv.png" into bracken_top_png
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv.svg" into bracken_top_svg

  script:
    """
    ${PYTHON3} \
    ${params.source_dir}/get-relative-abundance.py \
    --input-files ${reports} \
    --extension report \
    --classif-level ${level} \
    --output-dir . &&
    ${RSCRIPT} \
    ${params.source_dir}/plot-relative-abundance.Rscript \
    RES.${level}.rel_abundance.top_${params.max_species}.tsv \
    ${classified_reads} \
    ${level} \

    """
}

// compute diversity metrics

process diversity {

  cpus = 1
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/diversity", mode: "copy"

  input:
    tuple val(level), file(counts) from Bracken_counts

  output:
    file "RES.${level}.counts.diversity.tsv" into Diversity

  script:
    """
    ${RSCRIPT} \
    ${params.source_dir}/diversity-analysis.Rscript \
    ${counts} \
    . \
    ${level} \
    ${params.evenness} \
    ${params.dissimilarity} \

    """

}
