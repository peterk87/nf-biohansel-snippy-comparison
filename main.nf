#!/usr/bin/env nextflow

nextflow.preview.dsl=2

def helpMessage() {
  log.info"""
  ==================================================================
  ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
  ==================================================================

  Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

  Usage:
   The typical command for running the pipeline is as follows:
   
   nextflow run $workflow.manifest.name \\
     --outdir $params.outdir \\
     --schemesdir $params.schemesdir \\
     --n_genomes $params.n_genomes \\
     -work workdir \\
     -profile $workflow.profile

  Options:
    --outdir         Output directory (default: $params.outdir)
    --schemesdir     Directory with subtyping schemes and accessions to benchmark with biohansel (default: $params.schemesdir)
    --n_genomes      Number of SRA genomes to download and analyze per scheme (default: $params.n_genomes)
    --mincovs        List of minimum coverage values to test biohansel and snippy with delimited by comma (default: $params.mincovs)
  Other options:
    -w/--work-dir    The temporary directory where intermediate data will be saved (default: $workflow.workDir)
    -profile         Configuration profile to use. [singularity, conda, slurm] (default: $workflow.profile)
  Cluster options:
    -profile         Only "-profile slurm" is accepted
    --slurm_queue    Name of SLURM queue to submit jobs to (e.g. "HighPriority").
    --slurm_queue_size    SLURM queue size (default: $params.slurm_queue_size).
  """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  exit 1, "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

mincovs = []
if (params.mincovs instanceof String) {
  params.mincovs.split(',').each({
    try {
      mincovs << Integer.parseInt(it)
    } catch(Exception ex) {
      log.error "Error parsing mincovs integers: ${ex}"
    }
  })
} else if (params.mincovs instanceof Integer) {
  mincovs << params.mincovs
} else {
  exit 1, "Unexpected value '${params.mincovs}' for '--mincovs' of type ${params.mincovs.getClass()}"
}


if (mincovs.size() == 0) {
  exit 1, "You must specify at least one integer thread value in '--mincovs'"
}

// Header log info
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Schemes Directory'] = params.schemesdir
summary['# SRA genomes'] = params.n_genomes
summary['Min. coverages'] = mincovs
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "========================================="

outdir = params.outdir

baseSchemesDir = file(params.schemesdir)
if (!baseSchemesDir.exists()) {
  log.error "The schemes base directory does not exist at '$baseSchemesDir'!"
  exit 1
}

schemes = []

baseSchemesDir.eachDir { dir ->
  scheme = [dir.getName(), null, null, null]
  for (def f in dir.listFiles()) {
    switch(f.getName()) {
      case 'scheme.fasta':
        scheme[1] = f
        break
      case 'accessions':
        scheme[2] = f
        break
      case 'ref.gb':
        scheme[3] = f
        break
    }
  }
  schemes << scheme
}

if (schemes.size() == 0) {
  exit 1, "No biohansel/snippy analysis directories found in '${baseSchemesDir}'! Please specify a directory with one or more biohansel/snippy analysis directories."
}


process FASTERQ_DUMP {
  tag "$accession"
  publishDir "$outdir/fastqs/$scheme/$accession", mode: 'symlink', pattern: "*.fastq.gz"

  input:
    tuple val(scheme),
          path(scheme_fasta),
          val(accession),
          path(ref_genbank)
  output:
    tuple val(scheme),
          path(scheme_fasta),
          val(accession),
          path(reads1),
          path(reads2),
          path(ref_genbank)

  script:
  fq_1 = "${accession}_1.fastq"
  fq_2 = "${accession}_2.fastq"
  reads1 = "${fq_1}.gz"
  reads2 = "${fq_2}.gz"
  """
  fasterq-dump $accession -e ${task.cpus} -o $accession -S
  clumpify.sh -Xmx16g in=$fq_1 in2=$fq_2 out=$reads1 out2=$reads2 deleteinput=t
  """
}


process BIOHANSEL {
  tag "$scheme|$accession"
  publishDir "$outdir/biohansel/$accession", mode: 'copy', pattern: "*.tsv"

  input:
    val mincov
    tuple val(scheme),
          path(scheme_fasta),
          val(accession),
          path(reads1),
          path(reads2),
          path(ref_genbank)
  output:
    tuple val(scheme),
          val(accession),
          path(scheme_fasta),
          path(detailed_report),
          path(summary_report), emit: 'results'
    tuple val(scheme),
          val(1),
          val(1),
          val('biohansel_single'),
          path('.command.trace'),
          path(reads1),
          path(reads2), emit: 'trace'
  script:
  detailed_report = "biohansel-detailed_report-${accession}.tsv"
  summary_report = "biohansel-summary_report-${accession}.tsv"
  """
  mkdir -p reads
  ln -s `realpath *.fastq.gz` reads/
  hansel \\
    -v \\
    -t ${task.cpus} \\
    --min-kmer-freq $mincov \\
    -s $scheme_fasta \\
    -D reads/ \\
    -o $summary_report \\
    -O $detailed_report
  """
}


process SNIPPY {
  tag "$accession"
  publishDir "$outdir/snippy", mode: 'symlink'

  input:
    val mincov
    tuple val(scheme),
          path(scheme_fasta),
          val(accession),
          path(reads1),
          path(reads2),
          path(ref_genbank)
  output:
    tuple val(scheme),
          val(accession),
          path("$accession"), emit: 'results'
    tuple val(scheme),
          val(1),
          val(task.cpus),
          val('snippy'),
          path('.command.trace'),
          path(reads1),
          path(reads2), emit: 'trace'

  script:
  """
  snippy --prefix $accession \\
    --outdir $accession \\
    --cpus ${task.cpus} \\
    --ram ${task.memory.toGiga()} \\
    --mincov $mincov \\
    --R1 $reads1 \\
    --R2 $reads2 \\
    --ref $ref_genbank \\
    --tmpdir ./
  """
}


process COMPARE {
  tag "$accession"
  publishDir "$outdir/compare", pattern: "*.csv", mode: 'copy'

  input:
    tuple val(scheme), 
          val(accession), 
          path(scheme_fasta),
          path(bh_results), 
          path(bh_summary),
          path('snippy')
  output:
    path(csv)

  script:
  csv = "${accession}-${scheme}.csv"
  """
  compare_snippy_biohansel.py \\
    --bam-file snippy/${accession}.bam \\
    --biohansel-results $bh_results \\
    --snippy-consensus-subs-fa snippy/${accession}.consensus.subs.fa \\
    --scheme-fasta $scheme_fasta \\
    --output-csv $csv
  """
}



process TRACE_TABLE {
  publishDir "$outdir/trace", pattern: "*.csv", mode: 'copy'
  input:
    file trace
  output:
    file trace_table_csv

  script:
  trace_table_csv = "trace.csv"
  """
  make_trace_table.py -t $trace -o $trace_table_csv 
  """
}


workflow {


  Channel.from(schemes)
    .splitText(limit: params.n_genomes, elem: 2)
    .map { item ->
      item[2] = item[2].replaceAll("\\s", "")
      item
    }
    .filter { it[2] != '' }
    .dump(tag: "ch_accessions")
    .set { ch_accessions }

  mincov = Channel.value(mincovs.min())

  ch_accessions | FASTERQ_DUMP
  SNIPPY(mincov, FASTERQ_DUMP.out)
  BIOHANSEL(mincov, FASTERQ_DUMP.out)

  BIOHANSEL.out.results.join(SNIPPY.out.results, by: [0,1]) \
    | COMPARE

  BIOHANSEL.out.trace.mix(SNIPPY.out.trace)
    .collectFile() { scheme, samples, threads, type, trace_file, reads1, reads2 ->
      size_bytes = 0
      if (reads1 instanceof ArrayList) {
        reads1.collect( { size_bytes += file(it).size() } )
        reads2.collect( { size_bytes += file(it).size() } )
      } else {
        size_bytes = file(reads1).size() + file(reads2).size()
      }
      ['trace.txt', 
       """
       ${trace_file.text}
       scheme=${scheme}
       samples=${samples}
       threads=${threads}
       type=${type}
       size_bytes=${size_bytes}
       @@@
       """.stripIndent()]
    }
    .set { ch_trace }

    ch_trace | TRACE_TABLE
}

workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}


