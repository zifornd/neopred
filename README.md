[![GitHub Actions CI
Status](https://github.com/zifornd/neopred/actions/workflows/ci.yml/badge.svg)](https://github.com/zifornd/neopred/actions/workflows/ci.yml)
[![GitHub Actions Linting
Status](https://github.com/zifornd/neopred/actions/workflows/linting.yml/badge.svg)](https://github.com/zifornd/neopred/actions/workflows/linting.yml)[![Cite
with
Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with
conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with
docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with
singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera
Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/zifornd/neopred)

## Introduction

**neopred** is a bioinformatics pipeline that performs integrative
computational analysis of tumor immunity using bulk RNA-sequencing
(RNA-seq) data.

1.  Read QC
    ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2.  Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3.  Alignment ([`STAR`](https://github.com/alexdobin/STAR))
4.  Sort and index alignments
    ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5.  Quality Control
    ([`RSeQC`](https://github.com/MonashBioinformaticsPlatform/RSeQC))
6.  Gene Quantification
    ([`Salmon`](https://combine-lab.github.io/salmon/))
7.  Batch Removal
    ([`Limma`](https://www.bioconductor.org/packages/release/bioc/html/limma.html))
8.  Duplicate read marking
    ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
9.  Variant calling ([`GATK`](https://github.com/broadinstitute/gatk))
10. Variant Annotation
11. Neoantigen detection
    ([`arcasHLA`](https://github.com/RabadanLab/arcasHLA))
12. Prioritizing identified neoantigens
    ([`pVACseq`](https://github.com/griffithlab/pVAC-Seq))

## Usage

> [!NOTE] If you are new to Nextflow and nf-core, please refer to [this
> page](https://nf-co.re/docs/usage/installation) on how to set-up
> Nextflow. Make sure to [test your
> setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

````{=html}
First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness,batch,design,sample_group
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,auto,b1,control,WT_REP1_1
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse,b2,treatment,WT_REP2_1
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse,b1,control,SRR6357073_a,RAP1_UNINDUCED_REP1_1
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Rows with the same sample identifier are considered technical replicates and merged automatically. The strandedness refers to the library preparation and will be automatically inferred if set to auto.The batch refers to a parameter that evaluates if your samples have a batch effect. Design refers to the condition on which to do comparsion. The sample_group refers to condition used in neoanigen for removing duplicates.

-->
````

Now, you can run the pipeline using:

```bash
nextflow run neopred \
   -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING] Please provide pipeline parameters via the CLI or Nextflow
> `-params-file` option. Custom config files including those provided by
> the `-c` Nextflow option can be used to provide any configuration
> **_except for parameters_**; see
> [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage
documentation](https://nf-co.re/raredisease/usage) and the [parameter
documentation](https://nf-co.re/raredisease/parameters).

## Pipeline output {#pipeline-output}

To see the results of an example test run with a full size dataset refer
to the [results](https://nf-co.re/rnaseq/results) tab on the nf-core
website pipeline page. For more details about the output files and
reports, please refer to the [output
documentation](https://nf-co.re/raredisease/output).

## Credits

neopred was originally written by zifo.

We thank the following people for their extensive assistance in the
development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the
[contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can
be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by
the [nf-core](https://nf-co.re) community, reused here under the [MIT
license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics
> pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel,
> Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di
> Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi:
> [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
