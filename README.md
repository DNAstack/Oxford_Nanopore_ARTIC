# Oxford Nanopore SARS-CoV-2 data processing using the [ARTIC 
pipeline](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)

This repository provides a WDL wrapper for running the [Connor lab's implementation](https://github.com/connor-lab/ncov2019-artic-nf) of the [ARTIC pipeline](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) to process Oxford Nanopore single-ended SARS-CoV-2 sequencing data.


## Workflow inputs

An input template file with some defaults pre-defined can be found [here](./workflows/inputs.json).

| Input | Description |
|:-|:-|
| `accession` | Sample ID |
| `fastqs` | Array of single-ended FASTQ files |
| `primer_version` | The ARTIC primer version used to prepare the Nanopore library. Must be one of the schemes included in [the official ARTIC network repository](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019). Can be specified either as "ARTIC Vx" or as "Vx" |
| `min_length` | Minimum allowed read length for `artic guppyplex` [400] |
| `max_length` | Maximum allowed read length for `artic guppyplex` [700] |
| `min_reads` | Minimum number of reads output by `artic guppyplex` to pass QC. The workflow will exit with an error if this read count is not met. |
| `container_registry` | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |


## Workflow outputs

| Output | Description |
|:-|:-|
| `consensus_fa` | Genome assembly |
| `vcf`, `vcf_index` | Variant calls and index |
| `bam` | Reads aligned to the SARS-CoV-2 reference genome |
| `summary` | Pipeline metrics and plots |


## Containers

Docker image definitions can be found in our [bioinformatics-public-docker-images](https://github.com/DNAstack/bioinformatics-public-docker-images) repo.

The pipeline will always be pegged to a specific [ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) commit hash to avoid breaking due to updates to the underlying pipeline. The pipeline is periodically updated to the most recent version of the ncov2019-artic-nf pipeline.

All containers are publicly hosted in [DNAstack's container registry](https://hub.docker.com/u/dnastack).
