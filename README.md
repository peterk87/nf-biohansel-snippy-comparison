# nf-biohansel-snippy-comparison

[![Build Status](https://dev.azure.com/peterkruczkiewicz0831/nf-biohansel-snippy-comparison/_apis/build/status/peterk87.nf-biohansel-snippy-comparison?branchName=master)](https://dev.azure.com/peterkruczkiewicz0831/nf-biohansel-snippy-comparison/_build/latest?definitionId=2&branchName=master)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3444)


Nextflow workflow for comparing results from [biohansel](https://github.com/phac-nml/biohansel) and [Snippy](https://github.com/tseemann/snippy/) with NCBI SRA *Salmonella* Typhi genomes.


## Pre-reqs

- Install [Nextflow](https://www.nextflow.io/)
- Install [Conda](https://docs.conda.io/en/latest/miniconda.html)
- One or more directories each with the following files (see `schemes/typhi_v1.2` for an example)
 - `accessions` - List of SRA run accessions (e.g. `SRR8820085`) in a file (one accession per line)
 - `scheme.fasta` - biohansel scheme definition file
 - `ref.gb` - Genbank format reference genome
 
Input scheme directory included with this repo:

```
schemes
└── typhi_v1.2
    ├── accessions
    ├── ref.gb
    └── scheme.fasta
```

## Usage

Show help message:

```
nextflow run peterk87/nf-biohansel-snippy-comparison --help
```

Should see something like:

```
N E X T F L O W  ~  version 19.09.0-edge
Launching `main.nf` [wise_rubens] - revision: c9c6e63a09
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
==================================================================
peterk87/nf-biohansel-snippy-comparison  ~  version 1.0.0
==================================================================

Git info: null - null [null]

Usage:
 The typical command for running the pipeline is as follows:

 nextflow run peterk87/nf-biohansel-snippy-comparison \
   --outdir results \
   --schemesdir schemes \
   --n_genomes 5 \
   -work workdir \
   -profile standard

Options:
  --outdir         Output directory (default: results)
  --schemesdir     Directory with subtyping schemes and accessions to benchmark with biohansel (default: schemes)
  --n_genomes      Number of SRA genomes to download and analyze per scheme (default: 5)
  --mincovs        List of minimum coverage values to test biohansel and snippy with delimited by comma (default: 3,6,8)
Other options:
  -w/--work-dir    The temporary directory where intermediate data will be saved (default: ./work)
  -profile         Configuration profile to use. [singularity, conda, slurm] (default: standard)
Cluster options:
  -profile         Only "-profile slurm" is accepted
  --slurm_queue    Name of SLURM queue to submit jobs to (e.g. "HighPriority").
  --slurm_queue_size    SLURM queue size (default: 100).
```


Run test profile creating Conda environment:

```
nextflow run peterk87/nf-biohansel-snippy-comparison -profile test,conda
```

Run included benchmark dataset with Singularity and default parameters (i.e. first 5 genomes):

```
# clone/download this repo so that the scheme included with this repo can be run with the workflow
git clone https://github.com/peterk87/nf-biohansel-snippy-comparison.git
nextflow run peterk87/nf-biohansel-snippy-comparison -profile singularity --schemesdir nf-biohansel-snippy-comparison/schemes
```

Run above on a cluster with SLURM:

```
git clone https://github.com/peterk87/nf-biohansel-snippy-comparison.git
nextflow run peterk87/nf-biohansel-snippy-comparison -profile singularity,slurm --slurm_queue <QueueName> --schemesdir nf-biohansel-snippy-comparison/schemes
```

## Pipeline run information

Within your output directory (e.g. `results/`), you should find a `pipeline_info` directory with runtime information about your analysis including trace information (see https://www.nextflow.io/docs/latest/tracing.html for more info about these output files)

