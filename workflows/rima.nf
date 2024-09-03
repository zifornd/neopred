/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rima_pipeline'
include { getGenomeAttribute      } from '../subworkflows/local/utils_nfcore_rima_pipeline'
include { PREPARE_GENOME          } from '../subworkflows/local/prepare_genome'
include { PREPROCESS_STAR         } from '../subworkflows/local/preprocess_star'
include { QUANTIFY_SALMON         } from '../subworkflows/local/quantify_salmon'
include { PRE_VARIANTCALLING         } from '../subworkflows/local/pre_variantcalling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
//params.fasta = getGenomeAttribute('fasta')
//params.gtf = getGenomeAttribute('gtf')

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

workflow RIMA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input


    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Stage dummy file to be used as an optional input where required
    ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUBWORKFLOW: Prepare Genome
    //

    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.transcript_fasta
        )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Preprocess STAR
    //


    PREPROCESS_STAR (
        ch_samplesheet,
        PREPARE_GENOME.out.star_index.map { [ [:], it ] },
        PREPARE_GENOME.out.gtf.map { [ [:], it ] },
        params.star_ignore_sjdbgtf,
        '',
        params.seq_center ?: '',
        is_aws_igenome,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] }
    )

    ch_transcriptome_bam = PREPROCESS_STAR.out.bam_transcript
    ch_sorted_bam    = PREPROCESS_STAR.out.bam_sort
    ch_bam_bai       = PREPROCESS_STAR.out.bam_bai
    samtools_stats   = PREPROCESS_STAR.out.stats
    star_metrics     = PREPROCESS_STAR.out.metrics
    ch_versions      = ch_versions.mix(PREPROCESS_STAR.out.versions)

    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_STAR.out.log_final.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS_STAR.out.stats.collect{it[1]})

    QUANTIFY_SALMON (
        ch_transcriptome_bam,
        ch_dummy_file,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        true,
        params.salmon_quant_libtype ?: ''
    )

    ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)

    PRE_VARIANTCALLING(
        ch_bam_bai,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] },
        PREPARE_GENOME.out.fasta_fai.map { [ [:], it ] },
        [ [], 0 ]
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )



    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
