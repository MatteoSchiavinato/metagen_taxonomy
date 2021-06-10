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

  --ccs_dir               Directory containing FASTA files with ccs reads
  --extension             Extension of the FASTA files (fa)

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

Channel
  .fromPath("${params.ccs_dir}/*.${params.extension}")
  .map{ it -> [ it.simpleName, it] }
  .set{ Ccs }


// kraken2

process kraken2 {

  cpus = params.threads
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/kraken2", mode: "copy"

  input:
    tuple val(sample_id), file(fasta) from Ccs

  output:
    path "${sample_id}" into Kraken2_out
    tuple val(sample_id), file("${sample_id}/${sample_id}.report") into Counts

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

// bracken

process bracken {

  cpus = params.threads
  executor = "local"
  publishDir "${params.output_dir}/taxonomy/bracken", mode: "copy"

  input:
    tuple val(sample_id), file(report) from Counts

  output:
    path "${sample_id}" into Bracken_out
    tuple val("S"), path("S") into Bracken_S
    tuple val("G"), path("G") into Bracken_G
    tuple val("F"), path("F") into Bracken_F
    tuple val("C"), path("C") into Bracken_C
    tuple val("O"), path("O") into Bracken_O
    tuple val("P"), path("P") into Bracken_P
    tuple val("D"), path("D") into Bracken_D

  script:
    """
    if [ ! -d ${sample_id} ]; then mkdir ${sample_id}; fi &&
    unset X &&
    declare -a X=(S G F O C P D) &&
    for LEVEL in \${X[@]}
    do
      if [ ! -d \${LEVEL} ]; then mkdir ${sample_id}/\${LEVEL}; fi &&
      { ${BRACKEN} \
      -d ${params.kraken_db} \
      -i ${report} \
      -o ${sample_id}/\${LEVEL}/${sample_id}.\${LEVEL}.bracken \
      -w ${sample_id}/\${LEVEL}/${sample_id}.\${LEVEL}.bracken.report \
      -r ${params.kraken_db_read_len} \
      -l \${LEVEL} \
      -t ${params.min_counts}; } \
      &> ${sample_id}/\${LEVEL}/${sample_id}.bracken.log
    done \

    """
}


// combine channels

Bracken_S
  .mix(Bracken_G, Bracken_F, Bracken_C, Bracken_O, Bracken_P, Bracken_D)
  .set{ Bracken_reports }


// compute relative abundance

process plot_relative_abundance {

  cpus = 1
  executor = "local"
  publishDir "${params.output_dir}/taxonomy", mode: "copy"

  input:
    tuple val(level), path(level_reports) from Bracken_reports

  output:
    file "RES.${level}.rel_abundance.tsv" into bracken_frac
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv" into bracken_top
    file "RES.${level}.counts.tsv" into bracken_counts
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv.png" into bracken_top_png
    file "RES.${level}.rel_abundance.top_${params.max_species}.tsv.svg" into bracken_top_svg

  script:
    """
    ${PYTHON3} \
    ${params.source_dir}/get-relative-abundance.py \
    --input-dir ${level_reports} \
    --extension report \
    --classif-level ${level} \
    --output-dir . &&
    ${RSCRIPT} \
    ${params.source_dir}/plot-relative-abundance.Rscript \
    RES.${level}.rel_abundance.top_${params.max_species}.tsv \
    ${level} \

    """
}

// compute richness and diversity

// plot richness and diversity
