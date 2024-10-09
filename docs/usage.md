# nf-core/nf/rima: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/nf/rima/usage](https://nf-co.re/nf/rima/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## RUNNING THE PIPELINE

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-rima/ --input samplesheet.csv -profile docker --outdir <OUTDIR>
```
This will launch the pipeline with the `docker` configuration profile and save the outputs within respective folder inside `<OUTDIR>` folder. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                             # Directory containing the nextflow working files
<OUTDIR>                         # Finished results in specified location (defined with --outdir)
.nextflow_log                    # Hidden log file from Nextflow
```

If the user intends to run the pipeline using `conda` profile, the command to initiate the pipeline run is as follows:

```bash
nextflow run nf-rima/ --input samplesheet.csv -profile conda --outdir <OUTDIR>
```
This will launch the pipeline with the `conda` configuration profile and save the outputs within respective folder inside `<OUTDIR>` folder. The pipeline outputs will be created in the working directory and the output files generated will be as same as the run using the `-profile docker` option.

Similarly, apart from `docker` and `conda`, nf-rima pipeline is configured to run using different profiles namely, `mamba`, `singularity`, `podman`, `shifter`, `charliecloud`, `wave`, `apptainer`, `gitpod`, `test` and `test_full`. The user can run the pipeline utilizing different configured profiles by replacing `conda` in the above command with any of the above mentioned profile name.

To run pipeline using small test dataset provided in the pipeline repo with conda, use the command below:

```bash
nextflow run nf-rima/ -profile test,conda --outdir <OUTDIR>
```

To run pipeline using full complete test dataset provided in the pipeline repo with conda, use the command below:

```bash
nextflow run nf-rima/ -profile test_full,conda --outdir <OUTDIR>
```
Instead of the `conda` profile, to run the pipeline using docker profile for the small test dataset or full completed test dataset, the user will be required to replace `conda` with `docker` in the above command. 

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. The file supplied using `-params-file` option should be in `yaml` or `JSON` format. Moreover, pipeline specific parameters can be given via command line as `--<pipeline-parameter>` overwriting anyother values for the specific parameters passed via different configuration files. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration specific to compute environment except for parameters; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Samplesheet input

Create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```
Note, that the above input flag/parameter uses _double-hyphen_(`--`) as it is a pipeline parameter and not a nextflow parameter. The above input flag can be omitted if the full path to the file is mentioned in the nextflow.config and test.config files. The samplesheet has to be a comma-seperated file with at least 7 columns, and a header row as shown in the example below.

```console
sample,fastq_1,fastq_2,strandedness,design,batch,sample_group
SRR8281226,/home/saddam/samples/fastq/SRR8281226_chr6_sort_fixmate_dedup_R1.fq.gz,/home/saddam/samples/fastq/SRR8281226_chr6_sort_fixmate_dedup_R2.fq.gz,auto,control,b1,SRR8281226_1
SRR8281243,/home/saddam/samples/fastq/SRR8281243_chr6_sort_fixmate_dedup_R1.fq.gz,/home/saddam/samples/fastq/SRR8281243_chr6_sort_fixmate_dedup_R2.fq.gz,reverse,treatment,b2,SRR8281243_1
```

It is recommended to use gzipped fastq file with the extension `.fastq.gz` or `.fq.gz`. The `strandedness` column must be filled with small letters i.e. reverse and not Reverse/REVERSE. The samplesheet can have as many columns as you desire, however, there is a strict requirement for at least 7 columns to match those defined in the table below.

| Column              | Description                                                                                                                                                                            |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`            | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`           | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`           | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandedness`      | Sample strand-specificity. Must be one of `unstranded`, `forward` or `reverse`.                                                                                                        |
| `design`         | The name of the condition a sample belongs to (e.g. 'control', or 'treatment') - these labels will be used for downstream analysis.                                                    |
| `batch`        | The column indicates a group name indicating the samples belonging to the group name were processed together whereas a different group name indicate the sample was processed differently along with other set of samples.                                                                     |
| `sample_group` | The sample group name required for the downstream analysis .                                                       |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Reference Genome
The user may use a reference genome configured in the pipeline using iGenomes. To do that, the user may use `--genome <Reference-build>`. For example, the <Reference-build> could be either GRCh38 or GRCh37. if the user do not want to use `--genome` option, the user **must** provide to the pipeline a FASTA and a GTF filepath, without which the pipeline may will fail to run. The fasta and gtf filepath can be passed to the pipeline through config files such as nextflow.config or test_full.config file. Also the user can pass the filepath via commandline using pipeline parameters `--fasta </path/to/reference.fa.gz>` and `--gtf </path/to/reference.gtf.gz>`. For example, the user could save these reference file locally and provide the local file path to the pipeline parameters. The indices file required for STAR alignment is also created by the pipeline itself using the user provided reference fasta file. It should be noted, however, that the sequence files and indices of many common species are available for download from [AWS iGenomes](https://nf-co.re/usage/reference_genomes). 

# Preprocessing of FastQ reads

## FastQC
The tool [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a general QC standard for any sequencing workflow and nf-rima pipeline runs this at the beginning of the pipeline to generate general quality metrics of the sequenced reads. For example, FASTQC provides numerous quality control scores such as quality score distribution across reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For more details on parameters, please refer to [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

## Alignment 

The pipeline uses [STAR](https://github.com/alexdobin/STAR) for the **STAR_ALIGN** process to map raw FastQ reads to a reference genome and to project the alignments onto the transcriptome. STAR, by default, is configured in this pipeline to accept FASTQ files in gzipped format. The aligner tool requires fasta index file either to be supplied by the user via config file or fasta index file can be generated by the pipeline via the **GUNZIP_FASTA** process before the aligning the FastQ reads of a sample. To generate the index file, the user would require to provide the path to a reference annotation (GTF) file via the config file or using the parameter `--gtf </path/to/reference.gtf.gz>` containing the path to the local reference GTF file. If the index file are not provided by the user, **GUNZIP_FASTA** process is configured to run only once and designed to reuse the generated fasta index file for the samples requiring the subsequent alignment process.The alignment process use `--seq_center`, `--star_ignore_sjdbgtf`, and other optional alignment parameters, which the user can configure while running the star aligner. After the individual alignment runs of each sample, the pipeline via **STAR_METRICS** process will merge the alignment reports from all samples to generate a combined report. 

## Alignment Statistics
After alignment processing, the BAM files generated for each sample are processed using [samtools](https://github.com/samtools/samtools). Sorting of BAM files followed by indexing of BAM files and alignment statistics are generated for each of the BAM file.

## Alignment Quality Check
The pipeline use [RseQC](https://github.com/MonashBioinformaticsPlatform/RSeQC) to check the quality of read alignments and provide the principal measurements of RNA-seq data quality. To reduce the run time, the pipeline performs BAM downsampling and optionally (but recommended) if the path of bed file containing only the housekeeping gene regions is supplied via the config file or via command line using the pipeline parameter `--hk_bed </path/to/housekeeping.bed>`, then the housekeeping genes alignment alone are taken for the subsample. By default, the pipeline takes the bed file uploaded to the pipeline repositary, and without that parameter, pipeline can still run but at the cost of high run time required for downtream QC steps generating median Transcript Integrity Number (TIN) score, read distribution, and genebody coverage. The pipeline generates the quality metrics for an individual sample as well as combined quality metrics for the sets of samples used for running the workflow.

## Quantification
[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) tool is used to quantifiy the transcript level abundances from the transcriptome aligned BAM files. [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) set to run in alignment-based mode use the BAM file generated from the previous star process module and reference transcript fasta file as its input. Other optional additional arguments can be passed as parameter `--extra_salmon_quant_args` via the config file. Following [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) quantification process, the pipeline proceeds with **TX2GENE** and **TXIMPORT** process, to associate transcripts with geneID for gene-level summarization. In the **TX2GENE** process, the user need to set parameters `--gtf_extra_attributes` and `--gtf_group_features` appropriately for running other downstream analysis smoothly without any errors. 

## Batch effect removal
In order to avoid confounding the actual biological variation with the effects of the experimental design, the pipeline utilize **BATCH_REMOVAL** process to correct the batch effect between the samples used together for the analysis. The **BATCH_REMOVAL** module uses a Rscript [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) R package to remove any observed batch effect between the samples. In addition to the transcript quantification file from the salmon output, samplesheet CSV file, batch and design value need to be passed s input to the module.

## Somatic Variant calling - preprocessing of BAM
Variant calling workflow in the pipeline uses [GATK](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) best practices workflow. The sorted BAM file from the previous process will be used as a primary input, and in addition, preprocessing of BAM files occurs at each step of the variant calling workflow. In **PICARD_ADDORREPLACEREADGROUPS** process, read from the sorted BAM file are assigned to groups based on parameters supplied by the user. The read group library (RGPL) and read group platform unit (RGPU) value are fixed with default string value, whereas read group platform (RGPL) value as string can be changed by the user via the config (nextflow.config) file. By default, read group sample name (RGSM) value will be the respective sample id of that sample and using the optional parameter `--parse_read_group_info_file` the user can provide a text file containing samples and its associated RGSM value. The RGSM value for the specifc sample can be parsed and assigned to the respective sample. In **PICARD_MARKDUPLICATES** process, reads in a BAM file is located and duplicates from a single fragment of DNA are tagged. **PICARD_CREATESEQUENCEDICTIONARY** process is independent of the samples and generates a sequence ditcionary file from the reference genome FASTA file. The reference sequence dictionary file generated in the process will be used in the downstream process of the variant calling workflow. **GATK4_SPLITNCIGARREADS** split reads of the BAM file containing N CIGAR elements and utilizes marked duplicates BAM and FASTA file along with their indices. **GATK4_BASERECALIBRATOR** process, as a first pass in a two-stage process, recalibrates the quality score of each base by building an error model from the reads of the BAM file. The error model generation requires the user to supply, apart from the BAM and FASTA file, one or more known polymorphic sites (ExAc, gnomAD, or dbSNP resources) using the parameter `--known_indels </path/to/VCF-file>` as a VCF file. In addition, the pipeline requires an associated index file to be supplied via `--known_indels_index </path/to/VCF-Index-file>` parameter. **GATK4_APPLYBQSR**, as a second pass in a two-stage process, uses the recalibration table generated from BQSR step to recalibrate the base qualities of the input reads. 

## Somatic Variant Calling - variant identification & filtering
For calling somatic variants, **GATK4_MUTECT2** uses BAM file of the last process module as input along with FASTA, associated indices files, reference dictionary, and other optional supporting files. The user can optionally set minimum phred-scaled confidence threshold using the parameter `-emit-lod <user-entered-value>`. If another optional parameter `--F1R2_metrics` set to true, the process can generate an output *.f1r2.tar.gz file containing F1R2 counts, which is used for the downstream filtering step. If `--F1R2_metrics` set to false, the process will not generate the output file containing F1R2 counts. In addition, the variant identification process uses optional germline resource and panel of normal (PON) files for identifying variants via optional parameters `--germline_resource`, `--germline_resource_tbi`, `--pon` and `--pon_tbi` defined in the nextflow.config file. This pipeline is designed to give user an option to continue variant filter process or skip and directly proceed to variant annotation process. The user need to configure and set the optional parameter to true if the variant filtering steps are to be followed immediately after variant identification process. 

The BAM file from **GATK4_APPLYBQSR** step is passed to **GATK4_GETPILEUPSUMMARIES** process to generate metrics that can support or undermine contamination. For this purpose, in addition to BAM, FASTA, associated indices file and reference genome dictionary file, the user is also required to supply a common germline VCF file with popoulation allele frequencies using parameters `--pileup_vcf` and `--pileup_vcftbi` via the nextflow.config file. **GATK4_LEARNREADORIENTATIONMODEL** use f1r2 counts collected during mutect2 to learn the prior probability of read orientation artifacts and create the output *.tar.gz file. If `--F1R2_metrics` set to false, the mutect2 variant identification step will not generate *.f1r2.tar.gz file and the pipeline will skip this process step without creating the *.tar.gz file. In the next step, **GATK4_CALCULATECONTAMINATION** process, using pileup data from pileup summary process above, calculates the fraction of reads coming from cross-sample contamination. The user can set the optional parameter `--params.segmentation` to true, to generate an output table containing segmentation of tumor by minor allele fraction. The resulting contamination table, segmentation table, contamination estimate, and read orientation model ( *.f1r2.tar.gz file) are further used as optional input along with the required input VCF, FASTA, associated index and reference genome dict file are passed as input to **GATK4_FILTERMUTECTCALLS** process. The process filters SNV and indel from the VCF generated by mutect2 with without the optional files. **BCFTOOLS_VIEW** and **BCFTOOLS_INDEX** processes subset filtered vcf and create respective index file for variants that have passed the above filtering metric. The following processes **GATK4_SELECTVARIANTS** and **GATK4_COUNTVARIANTS** select only SNPs from the VCF and counts the number of SNP variants in the final VCF file. 

## Variant Annotation
The pipeline uses [VEP](https://github.com/Ensembl/ensembl-vep) tool for annotating filtered variants. The annotation workflow involves three process as part of the workflow. The **PVACTOOLS_INSTALLVEPPLUGIN** downloads the neccessary plugin required for downstream epitope prediction workflow. **ENSEMBLVEP_DOWNLOAD** following plugin download will download all the neccesary cache files required for the VEP tool annotation of the VCF files. The pipeline uses parameters `--vep_cache`, `--vep_genome_assembly`, `--vep_species` and `--vep_cache_version` that an user would be required to fill either via config file or via command line to run variant annotation pipeline without any error. By default, `--vep_cache` is set to null, `--vep_genome_assembly` is set as GRCh38, `--vep_species` is set to homo_sampiens and `--vep_cache_version` to 112. For example, If the user intend to change the vep tool's cache version or any other parameter, the user can update the new value for that parameter through the config file or via command line. 

## HLA typing workflow
This workflow uses [arcasHLA](https://github.com/RabadanLab/arcasHLA) tool to map RNAseq reads against the GRCh38 human genome without alternate contigs using the STAR aligner with default parameters. **ARCASHLA_EXTRACT** process, which uses input as the sorted BAM files, extract chromosome 6 reads and related HLA sequences. **ARCASHLA_GENOTYPE** use the extracted paired-end fastq files from the above step to predict the most likely genotype. **ARCASHLA_MERGE** process merge all JSONS produced in the genotyping step into a single table. **ARCASHLA_CONVERT** process changes alleles in the TSV file from its input form to a specified grouped nomenclature. The user need to set the pipeline parameters `--gtf_extra_attributes` set to null and `--gtf_group_features` as 'gene_name' either via config file or command line to run the command without any error. **ARCASHLA_PLOT** generates the plot as a png file consisting the different HLA frequencies. 

## Epitope prediciton workflow
This epitope prediction workflow depends on all the above three workflows i.e. Batch removal, HLA typing, and Variant annotation outputs of the above workflows are used as inputs to the epitope prediction workflow. This workflow involves a series of steps which integrates tumor mutation and expression data required for pvacseq tool to predict neoantigens that can raise immune response. **TX2GENE_PVACSEQ**, **TXIMPORT** and **BATCH_REMOVAL** processes helps expression data to be in compatible format with pvacseq tool. **VATOOLS_REFTRANSCRIPTMISMATCHREPORTER** process identify the mismatch variants (i.e. between the actual wildtype amino acid sequence and the expected amino acid sequence) from the input VCF using 'soft' as `--filter_val` parameter value and along with **BCFTOOLS_VIEW** process filter out the mismatch variants to a separate file. **VATOOLS_VCFEXPRESSIONANNOTATOR** process add the gene expression information to the annotaed mismatched removed VCF. The addition of gene expression information is written to the GX format field of the VCF which is passed to pvacseq tool for the neoantigen prediction. **PVACTOOLS_PVACSEQ** makes prediction for the transcript on the annotated variants by using three different MHC Class I binding affinity prediction algorithm and the pipeline is tested by using the MHCflurry, NetMHCcons, and MHCnuggetsII Class I prediction algorithm. The different algorithms can be changed or modified by user using the `--callers` pipeline parameter either via the nextflow.config file or through the command line. In addition, `--neoantigen_epitope1_lengths` have been tested with pvacseq tool default epitope length value which can also be specified by the user through the config file or via the command line. 

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/nf/rima
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/nf/rima releases page](https://github.com/nf-core/nf/rima/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_nf/rima:nf/rima:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_nf/rima:nf/rima:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_nf/rima:nf/rima:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_nf/rima:nf/rima:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_nf/rima:nf/rima:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).


## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```