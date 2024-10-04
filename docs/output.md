# nf/rima: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

-   [FastQC](#fastqc) - Raw read QC
-   Alignment.
    -   [STAR](#star) - Fast spliced aware genome alignment.
-   Alignment post-processing.
    -   [SAMtools](#samtools) - Sort and index alignments.
-   Quantification.
    -   [Salmon](#salmon) - Transcriptome quantification.
-   Quality Control
    -   [RSeQC](#rseqc) - Quality check of read alignments.
    -   [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline.
-   Batch Effect Removal.
    -   [Limma](#batch_removal) - R package to correct batch effects using a two-way ANOVA approach .
    -   [PCA](#pca) - Principal components analysis to validate that a batch effect has been removed.
-   HLA Typing
    -   [arcasHLA](#arcasHLA) - Predict HLA types for both MHC Class I & Class II from the bulk RNA-seq data.
-   [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC {#fastqc}

<details>

<summary>Output files</summary>

-   `fastqc/`
    -   `*_fastqc.html`: FastQC report containing quality metrics.
    -   `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

## Alignment {#alignment}

### STAR {#star}

<details>

<summary>Output files</summary>

-   `star/`
    -   `*.Aligned.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the reference genome will be placed in this directory.
    -   `*.Aligned.toTranscriptome.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the transcriptome will be placed in this directory.
-   `star/log/`
    -   `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
    -   `*.Log.final.out`: STAR alignment report containing the mapping results summary.
    -   `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
-   `star/unmapped/`
    -   `*.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

</details>

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for splice aware mapping typical of RNA sequencing data. STAR stands for *S*pliced *T*ranscripts *A*lignment to a *R*eference, and has been shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. Using `--aligner star_salmon` is the default alignment and quantification option.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as *Uniquely mapped* and few *Unmapped* reads.

![MultiQC - STAR alignment scores plot](mqc_star.png)

## Alignment post-processing {#alignment-post-processing}

### SAMtools {#samtools}

<details>

<summary>Output files</summary>

-   `<ALIGNER>/`
    -   `<SAMPLE>.sorted.bam`: If `--save_align_intermeds` is specified the original coordinate sorted BAM file containing read alignments will be placed in this directory.
    -   `<SAMPLE>.sorted.bam.bai`: If `--save_align_intermeds` is specified the BAI index file for the original coordinate sorted BAM file will be placed in this directory.
    -   `<SAMPLE>.sorted.bam.csi`: If `--save_align_intermeds --bam_csi_index` is specified the CSI index file for the original coordinate sorted BAM file will be placed in this directory.
-   `<ALIGNER>/samtools_stats/`
    -   SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_mapped.png)

![MultiQC - SAMtools mapped reads per contig plot](images/mqc_samtools_idxstats.png)

## Quantification {#quantification}

### Salmon {#salmon}

<details>

<summary>Output files</summary>

-   `salmon/`
    -   `<SAMPLE>/`
        -   `aux_info/`: Auxiliary info e.g. versions and number of mapped reads.
        -   `cmd_info.json`: Information about the Salmon quantification command, version and options.
        -   `lib_format_counts.json`: Number of fragments assigned, unassigned and incompatible.
        -   `libParams/`: Contains the file `flenDist.txt` for the fragment length distribution.
        -   `logs/`: Contains the file `salmon_quant.log` giving a record of Salmon's quantification.
        -   `quant.sf`: Salmon *transcript*-level quantification of the sample, including feature length, effective length, TPM, and number of reads.
-   `*lib_format_counts.json`: Contains information about the library format that was inferred.
-   `*meta_info.json`: Meta information from Salmon quant such as version and options used.

</details>

[Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying the expression of transcripts using RNA-seq data. It uses a quasi-mapping-based approach to provide fast and accurate quantification.

![MultiQC - Salmon fragment length distribution plot](images/mqc_salmon.png)

## Quality Control {#quality-control}

### RSeQC {#rseqc}

[RSeQC]((http://rseqc.sourceforge.net/)) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. You can tweak the supported scripts you would like to run by adjusting the `--rseqc_modules` parameter which by default will run all of the following: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

#### BAM down-sampling

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/down_sample/`
    -   `*_downsampling.bam`: File containing sub-sample of the original coordinate sorted BAM file.

</details>

This script sub-sample sorted BAM files to be used by RseQC to assess alignment quality.

#### Down-sampling Housekeeping

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/down_sample/`
    -   `*_downsample_hk.bam`: File containing alignments of house keeping genes for the downampled BAM file.

</details>

`bedtools intersect` allows the identification of overlaps between genomic features. In this context, the tool extracts the alignment of housekeeping genes using the input `BED` file.

#### TIN

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/tin/`
    -   `*.summary.txt`: File containing TIN results summary.
    -   `*.tin.xls`: XLS file containing TIN results.

</details>

This script is designed to evaluate RNA integrity at the transcript level. TIN (transcript integrity number) is named in analogous to RIN (RNA integrity number). RIN (RNA integrity number) is the most widely used metric to evaluate RNA integrity at sample (or transcriptome) level. It is a very useful preventive measure to ensure good RNA quality and robust, reproducible RNA sequencing. This process isn't run by default - please see [this issue](https://github.com/nf-core/rnaseq/issues/769).

RSeQC documentation: [tin.py](http://rseqc.sourceforge.net/#tin-py)

#### Read distribution

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/read_distribution/`
    -   `*.read_distribution.txt`: File containing fraction of reads mapping to genome feature e.g. CDS exon, 5'UTR exon, 3' UTR exon, Intron, Intergenic regions etc.

</details>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA-seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample but may be expected for a total RNA preparation.

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)

![MultiQC - RSeQC read distribution plot](images/mqc_rseqc_readdistribution.png)

#### Junction saturation

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/junction_saturation/pdf/`
    -   `*.junctionSaturation_plot.pdf`: PDF file containing junction saturation plot.
-   `<ALIGNER>/rseqc/junction_saturation/rscript/`
    -   `*.junctionSaturation_plot.r`: R script used to generate pdf plot above.

</details>

This script shows the number of splice sites detected within the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield any more observations. A good sample should approach such a plateau of *Known junctions*, however, very deep sequencing is typically required to saturate all *Novel Junctions* in a sample.

RSeQC documentation: [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)

![MultiQC - RSeQC junction saturation plot](images/mqc_rseqc_junctionsaturation.png)

#### Gene Body Coverage

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/gene_body_coverage/pdf`
    -   `*.geneBodyCoverage.curves.pdf`: PDF file containing gene body coverage curves.
-   `<ALIGNER>/rseqc/gene_body_coverage/rscript`
    -   `*.geneBodyCoverage.r`: R script used to generate pdf plot above.

</details>

This script calculates the RNA-seq read coverage over the gene body. Only sorted and indexed `BAM` files can be given as input. Genes or transcripts of mRNA length less than 100 are skipped.

#### TIN Summary

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/tinsummary`
    -   `*.tin_score_summary.txt`: File containing tin score summary.

</details>

This script summarizes the TIN score ranges across the samples.

#### Read Distribution matrix

<details>

<summary>Output files</summary>

-   `<ALIGNER>/rseqc/read_distributionmatrix`
    -   `read_distrib.matrix.tab`: A tab file containing distribution matrix across all samples.

</details>

This script creates RseQC read distribution matrix for all samples together.

### MultiQC {#multiqc}

<details>

<summary>Output files</summary>

-   `multiqc/`
    -   `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    -   `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    -   `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Batch Effect Removal {#batch-effect-removal}

### Limma {#limma}

<details>

<summary>Output files</summary>

-   `batch_removal/`
    -   `*_tpm.genesymbol.csv`: A csv file containing gene list across different samples before batch removal.
    -   `*_tpm.genesymbol.batchremoved.csv`: A csv file containing gene list across different samples after batch removal.

</details>

Batch effects across samples are easily overlooked but worth considering for immunotherapy cohort studies. Batch effects are usually caused by unbalanced experimental design and confound the estimation of group differences. To avoid confounding actual biological variation with the effects of experimental design, [`limma` ](https://doi.org/10.1093/nar/gkv007)and [`ComBat`](https://academic.oup.com/biostatistics/article/8/1/118/252073?login=true) are common approaches to correct batch effects. `Limma` uses a two-way ANOVA approach. `ComBat` uses an empirical Bayes approach, which is critical for small batches to avoid over-correction. For large batches, both methods should be similar. The [`sva`](https://bioconductor.org/packages/release/bioc/html/sva.html) R package implements both `ComBat` and surrogate variable analysis (sva) for batch effect correction.

### PCA {#pca}

<details>

<summary>Output files</summary>

-   `batch_removal/`
    -   `*_pca_plot_before.pdf`: PDF file containing PCA plot before batch removal.
    -   `*_pca_plot_after.pdf`: PDF file containing PCA plot after batch removal.

</details>

Principal components analysis (PCA) or unsupervised clustering before and after batch effect removal is an excellent way to validate that a batch effect has been removed. To evaluate if your samples have a batch effect, RIMA will generate PCA plots of gene expression data before and after batch effect removal by `limma`.

## HLA Typing

### arcas-HLA

### Pipeline information {#pipeline-information}

<details>

<summary>Output files</summary>

-   `pipeline_info/`
    -   Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    -   Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    -   Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
    -   Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
