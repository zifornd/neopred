# neopred: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- Alignment.
  - [STAR](#star) - Fast spliced aware genome alignment.
- Alignment post-processing.
  - [SAMtools](#samtools) - Sort and index alignments.
- Quantification.
  - [Salmon](#salmon) - Transcriptome quantification.
- Quality Control
  - [RSeQC](#rseqc) - Quality check of read alignments.
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline.
- Batch Effect Removal.
  - [Limma](#batch_removal) - R package to correct batch effects using a two-way ANOVA approach.
  - [PCA](#pca) - Principal components analysis to validate that a batch effect has been removed.
- HLA Typing
  - [arcasHLA](#arcasHLA) - Predict HLA types for both MHC Class I & Class II from the bulk RNA-seq data.
- [Variant Calling](#variant-calling)
  - [picard](#picard)
  - [GATK](#gatk)
- [Variant Annotation](#variant-annotation).
  - [VEP](#vep) - Transcriptome quantification.
- [pVACseq](#pvacseq)
  - [VAtools](#picard)
  - [pVACseq](#pvacseq)
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC {#fastqc}

<details>

<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

## Alignment {#alignment}

### STAR {#star}

<details>

<summary>Output files</summary>

- `star/`
  - `*.Aligned.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the reference genome will be placed in this directory.
  - `*.Aligned.toTranscriptome.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the transcriptome will be placed in this directory.
- `star/log/`
  - `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  - `*.Log.final.out`: STAR alignment report containing the mapping results summary.
  - `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
- `star/unmapped/`
  - `*.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

</details>

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for splice aware mapping typical of RNA sequencing data. STAR stands for *S*pliced *T*ranscripts *A*lignment to a *R*eference, and has been shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. Using `--aligner star_salmon` is the default alignment and quantification option.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

![MultiQC - STAR alignment scores plot](mqc_star.png)

## Alignment post-processing {#alignment-post-processing}

### SAMtools {#samtools}

<details>

<summary>Output files</summary>

- `<ALIGNER>/`
  - `<SAMPLE>.sorted.bam`: If `--save_align_intermeds` is specified the original coordinate sorted BAM file containing read alignments will be placed in this directory.
  - `<SAMPLE>.sorted.bam.bai`: If `--save_align_intermeds` is specified the BAI index file for the original coordinate sorted BAM file will be placed in this directory.
  - `<SAMPLE>.sorted.bam.csi`: If `--save_align_intermeds --bam_csi_index` is specified the CSI index file for the original coordinate sorted BAM file will be placed in this directory.
- `<ALIGNER>/samtools_stats/`
  - SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_mapped.png)

![MultiQC - SAMtools mapped reads per contig plot](images/mqc_samtools_idxstats.png)

## Quantification {#quantification}

### Salmon {#salmon}

<details>

<summary>Output files</summary>

- `salmon/`
  - `salmon.merged.gene_counts.tsv`: Matrix of gene-level raw counts across all samples.
  - `salmon.gene_tpm.tsv`: Matrix of gene-level TPM values across all samples.
  - `salmon.gene_counts.rds`: RDS object that can be loaded in R that contains a [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) container with the TPM (`abundance`), estimated counts (`counts`) and transcript length (`length`) in the assays slot for genes.
  - `salmon.merged.gene_lengths.tsv`: Matrix of average within-sample transcript lengths for each gene across all samples.
  - `salmon.merged.gene_counts_scaled.tsv`: Matrix of gene-level library size-scaled estimated counts across all samples.
  - `salmon.merged.gene_counts_length_scaled.tsv`: Matrix of gene-level length-scaled estimated counts across all samples.
  - `salmon.merged.transcript_counts.tsv`: Matrix of isoform-level raw counts across all samples.
  - `salmon.merged.transcript_tpm.tsv`: Matrix of isoform-level TPM values across all samples.
  - `tx2gene.tsv`: Tab-delimited file containing gene to transcripts ids mappings.

An additional subset of files for Salmon:

</details>

<details markdown="1">
<summary>Output files</summary>

- `salmon/`
  - `<SAMPLE>/`
    - `aux_info/`: Auxiliary info e.g. versions and number of mapped reads.
    - `cmd_info.json`: Information about the Salmon quantification command, version and options.
    - `lib_format_counts.json`: Number of fragments assigned, unassigned and incompatible.
    - `libParams/`: Contains the file `flenDist.txt` for the fragment length distribution.
    - `logs/`: Contains the file `salmon_quant.log` giving a record of Salmon's quantification.
    - `quant.sf`: Salmon _transcript_-level quantification of the sample, including feature length, effective length, TPM, and number of reads.
- `*lib_format_counts.json`: Contains information about the library format that was inferred.
- `*meta_info.json`: Meta information from Salmon quant such as version and options used.

</details>

[Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying the expression of transcripts using RNA-seq data. It uses a quasi-mapping-based approach to provide fast and accurate quantification.

![MultiQC - Salmon fragment length distribution plot](images/mqc_salmon.png)

## Quality Control {#quality-control}

### RSeQC {#rseqc}

[RSeQC](<(http://rseqc.sourceforge.net/)>) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. You can tweak the supported scripts you would like to run by adjusting the `--rseqc_modules` parameter which by default will run all of the following: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

#### BAM down-sampling

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/down_sample/`
  - `*_downsampling.bam`: File containing sub-sample of the original coordinate sorted BAM file.

</details>

This script sub-sample sorted BAM files to be used by RseQC to assess alignment quality.

#### Down-sampling Housekeeping

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/down_sample/`
  - `*_downsample_hk.bam`: File containing alignments of house keeping genes for the downampled BAM file.

</details>

`bedtools intersect` allows the identification of overlaps between genomic features. In this context, the tool extracts the alignment of housekeeping genes using the input `BED` file.

#### TIN

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/tin/`
  - `*.summary.txt`: File containing TIN results summary.
  - `*.tin.xls`: XLS file containing TIN results.

</details>

This script is designed to evaluate RNA integrity at the transcript level. TIN (transcript integrity number) is named in analogous to RIN (RNA integrity number). RIN (RNA integrity number) is the most widely used metric to evaluate RNA integrity at sample (or transcriptome) level. It is a very useful preventive measure to ensure good RNA quality and robust, reproducible RNA sequencing. This process isn't run by default - please see [this issue](https://github.com/nf-core/rnaseq/issues/769).

RSeQC documentation: [tin.py](http://rseqc.sourceforge.net/#tin-py)

#### Read distribution

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_distribution/`
  - `*.read_distribution.txt`: File containing fraction of reads mapping to genome feature e.g. CDS exon, 5'UTR exon, 3' UTR exon, Intron, Intergenic regions etc.

</details>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA-seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample but may be expected for a total RNA preparation.

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)

![MultiQC - RSeQC read distribution plot](images/mqc_rseqc_readdistribution.png)

#### Junction saturation

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/junction_saturation/pdf/`
  - `*.junctionSaturation_plot.pdf`: PDF file containing junction saturation plot.
- `<ALIGNER>/rseqc/junction_saturation/rscript/`
  - `*.junctionSaturation_plot.r`: R script used to generate pdf plot above.

</details>

This script shows the number of splice sites detected within the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield any more observations. A good sample should approach such a plateau of _Known junctions_, however, very deep sequencing is typically required to saturate all _Novel Junctions_ in a sample.

RSeQC documentation: [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)

![MultiQC - RSeQC junction saturation plot](images/mqc_rseqc_junctionsaturation.png)

#### Gene Body Coverage

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/gene_body_coverage/pdf`
  - `*.geneBodyCoverage.curves.pdf`: PDF file containing gene body coverage curves.
- `<ALIGNER>/rseqc/gene_body_coverage/rscript`
  - `*.geneBodyCoverage.r`: R script used to generate pdf plot above.

</details>

This script calculates the RNA-seq read coverage over the gene body. Only sorted and indexed `BAM` files can be given as input. Genes or transcripts of mRNA length less than 100 are skipped.

#### TIN Summary

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/tinsummary`
  - `*.tin_score_summary.txt`: File containing tin score summary.

</details>

This script summarizes the TIN score ranges across the samples.

#### Read Distribution matrix

<details>

<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_distributionmatrix`
  - `read_distrib.matrix.tab`: A tab file containing distribution matrix across all samples.

</details>

This script creates RseQC read distribution matrix for all samples together.

### MultiQC {#multiqc}

<details>

<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Batch Effect Removal {#batch-effect-removal}

### Limma {#limma}

<details>

<summary>Output files</summary>

- `batch_removal/`
  - `*_tpm.genesymbol.csv`: A csv file containing gene list across different samples before batch removal.
  - `*_tpm.genesymbol.batchremoved.csv`: A csv file containing gene list across different samples after batch removal.

</details>

Batch effects across samples are easily overlooked but worth considering for immunotherapy cohort studies. Batch effects are usually caused by unbalanced experimental design and confound the estimation of group differences. To avoid confounding actual biological variation with the effects of experimental design, [`limma` ](https://doi.org/10.1093/nar/gkv007)and [`ComBat`](https://academic.oup.com/biostatistics/article/8/1/118/252073?login=true) are common approaches to correct batch effects. `Limma` uses a two-way ANOVA approach. `ComBat` uses an empirical Bayes approach, which is critical for small batches to avoid over-correction. For large batches, both methods should be similar. The [`sva`](https://bioconductor.org/packages/release/bioc/html/sva.html) R package implements both `ComBat` and surrogate variable analysis (sva) for batch effect correction.

### PCA {#pca}

<details>

<summary>Output files</summary>

- `batch_removal/`
  - `*_pca_plot_before.pdf`: PDF file containing PCA plot before batch removal.
  - `*_pca_plot_after.pdf`: PDF file containing PCA plot after batch removal.

</details>

Principal components analysis (PCA) or unsupervised clustering before and after batch effect removal is an excellent way to validate that a batch effect has been removed. To evaluate if your samples have a batch effect, neopred will generate PCA plots of gene expression data before and after batch effect removal by `limma`.

## HLA Typing

### arcasHLA

#### arcasHLA extract

<details>

<summary>Output files</summary>

- `arcashla/extract`
  - `<SAMPLE>.extracted.1.fq.gz,<SAMPLE>.extracted.2.fq.gz`: extracted chromosome 6 reads and related HLA sequences, By default, extract outputs paired FASTQ files.
  - `<SAMPLE>.log` : log file for run summary.

</details>

This module extracts reads mapped to chromosome 6 and any HLA decoys or chromosome 6 alternates.

#### arcasHLA genotype

<details>

<summary>Output files</summary>

- `arcashla/genotype`
  - `<SAMPLE>.alignment.p`: Contains alignment information for the reads, crucial for HLA typing.
  - `<SAMPLE>.genes.json` : Stores gene-level information from the analysis, detailing HLA genes.
  - `<SAMPLE>.genotype.json` : Holds the genotype results, indicating the HLA alleles identified.
  - `<SAMPLE>.genotype.log` : A log file recording the steps and status of the genotyping process.

</details>

This module genotypes HLA alleles from extracted reads (no partial alleles).

#### arcasHLA merge

<details>

<summary>Output files</summary>

- `arcashla/merge`
  - `genotypes.tsv`: : This file contains the combined HLA genotyping results from multiple samples in a tab-separated values format.

</details>

This module merge genotyping output for multiple samples into a single json file.

#### arcasHLA convert

<details>

<summary>Output files</summary>

- `arcashla/convert`
  - `genotypes.p-group.tsv`: This file contains HLA genotyping results converted to P-group nomenclature, making it easier to compare and standardize HLA data across different datasets.

</details>

This script sub-sample sorted BAM files to be used by RseQC to assess alignment quality.

## Variant Calling {#variant-calling}

### picard {#picard}

#### picard createsequencedictionary

<details>

<summary>Output files</summary>

- `picard/createsequencedictionary`
  - `<GENOME>.genome.dict`: : A file that contains a sequence dictionary for the reference genome.

</details>

[CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard): Creates a sequence dictionary for a reference sequence in FASTA format.

#### picard collectmultiplemetrics

<details>

<summary>Output files</summary>

- `picard/picard_metrics`
  - `<SAMPLE>.mLb.clN.CollectMultipleMetrics.alignment_summary_metrics` : Contains alignment summary statistics, detailing how well reads are aligned to the reference genome.
  - `<SAMPLE>.mLb.clN.CollectMultipleMetrics.base_distribution_by_cycle_metrics` : Provides base distribution information by cycle, showing the frequency of each base (A, T, C, G) at each sequencing cycle.
  - `<SAMPLE>.mLb.clN.CollectMultipleMetrics.insert_size_metrics` : Contains metrics about the insert sizes of paired-end reads, useful for assessing library preparation quality.
  - `<SAMPLE>.mLb.clN.CollectMultipleMetrics.quality_by_cycle_metrics` : Shows the quality scores of bases by sequencing cycle, indicating sequencing accuracy across cycles.
  - `<SAMPLE>.mLb.clN.CollectMultipleMetrics.quality_distribution_metrics` : Provides quality score distributions, showing the overall quality of the sequencing reads.
  - `<SAMPLE>.markdup.sorted.MarkDuplicates.metrics.txt` : Contains metrics from the MarkDuplicates tool, detailing the number of duplicate reads detected and removed, which helps in evaluating library complexity.

</details>

[CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard): Collects multiple classes of metrics from a SAM/BAM file in a single pass.

#### picard addorreplacereadgroups

<details>

<summary>Output files</summary>

- `picard/addorreplacereadgroups`
  - `<SAMPLE>.addorreplacereadgroups.bam`: This is the BAM file with updated or added read group information.
  - `<SAMPLE>.addorreplacereadgroups.bai` : The index file for the BAM.

</details>

[AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard): Adds or replaces read groups in a BAM file.

#### picard markduplicates

<details>

<summary>Output files</summary>

- `picard/markduplicates`
  - `<SAMPLE>.markdup.sorted.bam` : This BAM file contains alignments with duplicate reads marked.

</details>

[MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard): Identifies and tags duplicate reads in a BAM or SAM file.

#### samtools stats

<details>

<summary>Output files</summary>

- `picard/samtools_stats`
  - `<SAMPLE>.markdup.sorted.bam.flagstat` : Contains summary statistics about the marked duplicate BAM file.
  - `<SAMPLE>.readgroup.sorted.bam.flagstat` : Contains summary statistics about the read group sorted BAM file.

</details>

[samtools stats](http://www.htslib.org/doc/samtools-stats.html): Produces comprehensive statistics from a BAM file.

### GATK {#gatk}

#### gatk4 splitncigarreads

<details>

<summary>Output files</summary>

- `gatk4/splitncigarreads`
  - `<SAMPLE>.bam`: This BAM file contains the reads that have been split at introns, which is particularly useful for handling RNA-Seq data in downstream analyses, ensuring accurate variant calling and alignment.

</details>

[SplitNCigarReads](https://gatk.broadinstitute.org/hc/en-us/articles/360036450212-SplitNCigarReads): Splits reads with Ns in their CIGAR strings, typically for RNA-seq data.

#### gatk4 applybqsr

<details>

<summary>Output files</summary>

- `gatk4/applybqsr`
  - `<SAMPLE>.applybqsr.bam` : This BAM file contains the reads that have undergone Base Quality Score Recalibration (BQSR).

</details>

[ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR): Applies base quality score recalibration to a BAM or CRAM file.

#### gatk4 mutect2

<details>

<summary>Output files</summary>

- `gatk4/mutect2/<SAMPLE>`
  - `<SAMPLE>.f1r2.tar.gz`: Contains counts of F1 and R2 reads used to learn the read orientation model for filtering artifacts.
  - `<SAMPLE>.vcf.gz` : A gzipped VCF file with variant calls from Mutect2.
  - `<SAMPLE>.vcf.gz.tbi` : An index file for the gzipped VCF, enabling efficient access to the variant data.

</details>

[Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Calls somatic SNVs and indels via local assembly of haplotypes.

#### gatk4 getpileupsummaries

<details>

<summary>Output files</summary>

- `gatk4/getpileupsummaries`
  - `<SAMPLE>.pileups.table`: table that summarizes counts of reads supporting reference, alternate, and other alleles at given sites.

</details>

[GetPileupSummaries](https://gatk.broadinstitute.org/hc/en-us/articles/360037593451-GetPileupSummaries): Summarizes counts of reads supporting reference, alternate, and other alleles for given sites.

#### gatk4 learnreadorientationmodel

<details>

<summary>Output files</summary>

- `gatk4/learnreadorientationmodel`
  - `<SAMPLE>.f1r2.tar.gz`: file, which contains the learned artifact priors, These priors are then used by Mutect2 to filter out potential artifacts during variant calling.

</details>

[LearnReadOrientationModel](https://gatk.broadinstitute.org/hc/en-us/articles/360051305331-LearnReadOrientationModel): Learns the prior probability of read orientation artifacts from F1R2 counts.

#### gatk4 calculatecontamination

<details>

<summary>Output files</summary>

- `gatk4/calculatecontamination/<SAMPLE>`
  - `<SAMPLE>.contamination.table`: provides an estimate of the contamination fraction in the tumor sample.
  - `<SAMPLE>.segmentation.table` : this table is used to describe the segments of the genome that have been analyzed for contamination.

</details>

[CalculateContamination](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Estimates contamination in a sample using pileup summaries.

#### gatk4 filtermutectcalls

<details>

<summary>Output files</summary>

- `gatk4/filtermutectcalls`
  - `<SAMPLE>.filtered.vcf.gz`: containing only the filtered variants filtered using contamination and sementation table.

</details>

[FilterMutectCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Filters somatic variant calls made by Mutect2.

#### gatk4 selectvariants

<details>

<summary>Output files</summary>

- `gatk4/selectvariants`
  - `<SAMPLE>.selected.vcf.gz`: contains a subset of variants selected based on specified criteria, retaining the VCF format with the desired variants.

</details>

[SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Selects a subset of variants from a VCF file.

#### gatk4 countvariants

<details>

<summary>Output files</summary>

- `gatk4/countvariants`
  - `<SAMPLE>_counts.txt`: file have count of variant records in the specified VCF file.

</details>

[CountVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2): Counts the number of variants in a VCF file.

## Variant Annotation {#variant-annotation}

### VEP {#vep}

<details>

<summary>Output files</summary>

- `ensemblvep/vep`
  - `<SAMPLE>.annotated.vcf.gz` : VCF file provides detailed annotations of genetic variants, including their effects on genes, transcripts, proteins, and regulatory regions
  - `<SAMPLE>.annotated.vcf.gz_summary.html`: HTML file containing statistics of VEP output.

</details>

VEP: Variant Effect Predictor (VEP) annotates variants with information about their effects on genes and proteins.

## Epitope Prediction {#epitope-prediction}

### VAtools {#vatools}

<details>

<summary>Output files</summary>

- `vatools`
  - `<SAMPLE>.filter.vcf` : A file contain only the variants that match between the reference genome and the Ensembl reference transcript, ensuring accurate consequence annotations.
  - `<SAMPLE>.mutect2.somatic.base.snp.Somatic.hc.filter.vep.gx.vcf`: VCF file with added expression data from the specified tool, annotated in the INFO column.

</details>

[VCF Annotation Tools](https://github.com/griffithlab/VAtools) is a python package that includes several tools to annotate VCF files with data from other tools.

### pVACseq {#pvacseq}

<details>

<summary>Output files</summary>

- `pvacseq/MHC_Class_I`
  - `<SAMPLE>.tsv` : An intermediate file with variant, transcript, coverage, vaf, and expression information parsed from the input files.
  - `<SAMPLE>.tsv_<chunks>`: The above file but split into smaller chunks for easier processing with IEDB.
  - `<SAMPLE>.fasta` : A fasta file with mutant and wildtype peptide subsequences for all processable variant-transcript combinations.
  - `<SAMPLE>.all_epitopes.tsv`: A list of all predicted epitopes and their binding affinity scores, with additional variant information.
  - `<SAMPLE>.filtered.tsv` : The above file after applying all filters, with (optionally) cleavage site, stability predictions, and reference proteome similarity metrics added.
  - `<SAMPLE>.all_epitopes.aggregated.tsv`: An aggregated version of the all_epitopes.tsv file that gives information about the best epitope for each mutation in an easy-to-read format.
  - `ui.R, app.R, server.R, styling.R, anchor_and_helper_functions.R` : pVACview R Shiny application files.
  - `www` : Directory containing image files for pVACview. Not generated when running with elution algorithms only.

</details>

[pVAC-Seq](https://github.com/griffithlab/pVAC-Seq) a flexible, streamlined computational workflow for identification of personalized Variant Antigens by Cancer Sequencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq).

### Pipeline information {#pipeline-information}

<details>

<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
