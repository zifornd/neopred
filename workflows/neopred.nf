/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS : Consists of validation based plugins and modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS : Consists of a mix of local and nf-core subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_neopred_pipeline'
include { getGenomeAttribute      } from '../subworkflows/local/utils_nfcore_neopred_pipeline'
include { PREPARE_GENOME          } from '../subworkflows/local/prepare_genome'
include { PREPROCESS_STAR         } from '../subworkflows/local/preprocess_star'
include { RSEQC                   } from '../subworkflows/local/rseqc'
include { QUANTIFY_SALMON         } from '../subworkflows/local/quantify_salmon'
include { PRE_VARIANTCALLING      } from '../subworkflows/local/pre_variantcalling'
include { BATCH_REMOVAL_ANALYSIS  } from '../subworkflows/local/batch_removal_analysis'
include { HLA_TYPING              } from '../subworkflows/local/hla_typing'
include { VARIANT_CALLINGFILTERING } from '../subworkflows/local/gatk_varcall'
include { VARIANT_ANNOTATION      } from '../subworkflows/local/variant_annotation'
include { EPITOPE_PREDICTION      } from '../subworkflows/local/epitope_prediction'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

workflow NEOPRED {

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

    if ((!params.variant_calling) && (params.variant_filtering) && (!params.variant_annotation)) {

        error "ERROR: VARIANT FILTERING STEP REQUIRES VARIANT CALLING WORKFLOW AND VARIANT FILTERING ALONE CANNOT BE TRUNED ON \n" +
            "HINT: RUN VARIANT CALLING WORKFLOW FOLLOWED BY VARIANT FILTERING PROCESS AND VARIANT ANNOTATION WORKFLOW \n"

        //exit 1

        }

    if ((!params.variant_calling) && (params.variant_filtering) && (params.variant_annotation)) {

        error "ERROR: VARIANT FILTERING STEP AND VARIANT ANNOTATION REQUIRES RUNNING VARIANT CALLING WORKFLOW \n" +
            "HINT: RUN VARIANT CALLING WORKFLOW FOLLOWED BY VARIANT FILTERING PROCESS AND VARIANT ANNOTATION WORKFLOW \n"

        }

    if ((!params.variant_calling) && (!params.variant_filtering) && (params.variant_annotation)) {

        error "ERROR: VARIANT ANNOTATION REQUIRES VARIANT CALLING WORKFLOW TO BE TURNED ON \n" +
            "HINT: RUN VARIANT CALLING WORKFLOW FOLLOWED BY VARIANT ANNOTATION WORKFLOW \n"

        }

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

    //
    // SUBWORKFLOW: RSeQC
    //

    if (!params.hk_bed) {

        ch_hk_bed = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)

        RSEQC (
            ch_bam_bai,
            samtools_stats,
            params.gtf,
            ch_hk_bed
        )

        ch_tin_multiqc                = RSEQC.out.tin_txt.collect{it[1]}
        ch_tin_multiqc                = ch_tin_multiqc.mix(RSEQC.out.tin_summary.collect{it[1]})
        ch_junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript.collect{it[1]}
        ch_readdistribution_multiqc   = RSEQC.out.readdistribution_txt.collect{it[1]}
        ch_readdistribution_multiqc   = ch_readdistribution_multiqc.mix(RSEQC.out.readdistribution_matrix)
        ch_multiqc_files              = ch_multiqc_files.mix(ch_tin_multiqc,ch_junctionsaturation_multiqc,ch_readdistribution_multiqc)
        ch_down_bam_bai               = RSEQC.out.down_bam_bai
        //ch_hk_bam_bai                 = RSEQC.out.hk_bam_bai
        ch_versions                   = ch_versions.mix(RSEQC.out.versions)

    }

    else {

        ch_hk_bed = params.hk_bed

        RSEQC (
            ch_bam_bai,
            samtools_stats,
            params.gtf,
            ch_hk_bed
        )

        ch_tin_multiqc                = RSEQC.out.tin_txt.collect{it[1]}
        ch_tin_multiqc                = ch_tin_multiqc.mix(RSEQC.out.tin_summary.collect{it[1]})
        ch_junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript.collect{it[1]}
        ch_readdistribution_multiqc   = RSEQC.out.readdistribution_txt.collect{it[1]}
        ch_readdistribution_multiqc   = ch_readdistribution_multiqc.mix(RSEQC.out.readdistribution_matrix)
        ch_multiqc_files              = ch_multiqc_files.mix(ch_tin_multiqc,ch_junctionsaturation_multiqc,ch_readdistribution_multiqc)
        ch_down_bam_bai               = RSEQC.out.down_bam_bai
        //ch_hk_bam_bai               = RSEQC.out.hk_bam_bai
        ch_versions                   = ch_versions.mix(RSEQC.out.versions)
    }

    //
    // SUBWORKFLOW: Salmon Quantification
    //

    QUANTIFY_SALMON (
        ch_transcriptome_bam,
        ch_dummy_file,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        true,
        params.salmon_quant_libtype ?: ''
    )

    ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_SALMON.out.multiqc.collect{it[1]})

    //
    // SUBWORKFLOW: Batch removal and PCA
    //

    BATCH_REMOVAL_ANALYSIS (
        params.input,
        params.batch,
        params.design,
        QUANTIFY_SALMON.out.tpm_gene)

    ch_versions = ch_versions.mix(BATCH_REMOVAL_ANALYSIS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BATCH_REMOVAL_ANALYSIS.out.before_br_pca)
    ch_multiqc_files = ch_multiqc_files.mix(BATCH_REMOVAL_ANALYSIS.out.after_br_pca)

    if (params.arcasHLA) {
        //
        // SUBWORKFLOW: arcasHLA Typing
        //
        HLA_TYPING (
            params.input,
            ch_sorted_bam,
            BATCH_REMOVAL_ANALYSIS.out.tpm,
            params.batch,
            params.design,
            params.patient_id)

        ch_multiqc_files = ch_multiqc_files.mix(HLA_TYPING.out.hla_log.collect{it[1]})
        //ch_multiqc_files = ch_multiqc_files.mix(HLA_TYPING.out.hla_plot)
        ch_versions = ch_versions.mix(HLA_TYPING.out.versions)

    }

    if (params.variant_calling) {

        PRE_VARIANTCALLING(
            ch_sorted_bam,
            ch_bam_bai,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fasta_fai.map { [ [:], it ] },
            params.dbsnp,
            params.dbsnp_tbi
        )

        ch_bqsr_bam = PRE_VARIANTCALLING.out.bqsr_bam
        ch_bqsr_bai = PRE_VARIANTCALLING.out.bqsr_bai
        ch_dict     = PRE_VARIANTCALLING.out.dict
        ch_versions = ch_versions.mix(PRE_VARIANTCALLING.out.versions)

        //
        // Subworkflow: Variant Identification and Filtering using GATK
        //

        ch_fasta = PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        ch_fai   = PREPARE_GENOME.out.fasta_fai.map { [ [:], it ] }

        VARIANT_CALLINGFILTERING (
            ch_bqsr_bam,
            ch_bqsr_bai,
            ch_fasta,
            ch_fai,
            ch_dict,
            params.germline_resource,
            params.germline_resource_tbi,
            params.pon,
            params.pon_tbi,
            params.pileup_vcf,
            params.pileup_vcftbi
        )

        ch_variants         =  VARIANT_CALLINGFILTERING.out.selected_vcf
        //ch_f1r2             =  VARIANT_CALLINGFILTERING.out.f1r2
        ch_variants_tbi     =  VARIANT_CALLINGFILTERING.out.selected_tbi
        //ch_variants_stats   =  VARIANT_CALLINGFILTERING.out.variants_stats
        ch_versions         = ch_versions.mix( VARIANT_CALLINGFILTERING.out.versions)
}
    if ((params.variant_calling) && (params.variant_annotation)) {
            //
            // Subworkflow: Variant Annotation using VEP
            //

            // Download cache if needed
            // Assuming that if the cache is provided, the user has already downloaded it
            ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome_assembly}" ], params.vep_genome_assembly, params.vep_species, params.vep_cache_version ])
            //var=Channel.fromPath(params.raw_vcf)

            VARIANT_ANNOTATION (
                ch_variants,
                ch_variants_tbi,
                ensemblvep_info,
                ch_fasta,
                params.vep_genome_assembly,
                params.vep_species,
                params.vep_cache_version
            )

            ch_annot_vcf = VARIANT_ANNOTATION.out.results
            ch_versions         = ch_versions.mix( VARIANT_ANNOTATION.out.versions)
        }

    if ((!params.variant_calling) && (!params.variant_annotation) && (!params.arcasHLA)) {

        log.warn "EPITOPE PREDICTION PROCESS REQUIRES ARCASHLA MODULE, VARIANT CALLING AND VARIANT ANNOTATION MODULES TO BE TURNED ON \n"
        log.info "IF YOU NEED TO RUN EPITOPE PREDICTION MODULE, PLEASE TURN ON ARCASHLA MODULE, VARIANT CALLING AND VARIANT ANNOTATION MODULES \n"

    }

    if ((params.variant_calling) && (params.variant_annotation) && (!params.arcasHLA)) {

        log.warn "WARNING: EPITOPE PREDICTION PROCESS REQUIRES ARCASHLA MODULE TO BE TURNED ON \n"
        log.info "SUGGESTION: PLEASE TURN ON ARCASHLA MODULE MODULE TOO \n"

    }

    if ((params.variant_calling) && (!params.variant_annotation) && (params.arcasHLA)) {

        log.warn "WARNING: EPITOPE PREDICTION PROCESS REQUIRES VARIANT ANNOTATION MODULE TO BE TURNED ON \n"
        log.info "SUGGESTION: PLEASE TURN ON VARIANT ANNOTATION MODULE TOO"

    }

    if ((params.variant_calling) && (!params.variant_annotation) && (!params.arcasHLA)) {

        log.warn "WARNING: EPITOPE PREDICTION PROCESS REQUIRES ARCASHLA MODULE AND VARIANT ANNOTATION MODULES TO BE TURNED ON \n"
        log.info "SUGGESTION: PLEASE TURN ON ARCASHLA MODULE FOR RUNNING EPITOPE PREDICTION WORKFLOW \n"

    }

    if ((params.variant_calling) && (params.variant_annotation) && (params.arcasHLA)) {

    //
    // Subworkflow: Epitope prediction using pVACseq
    //

        pvacseq_geno =  HLA_TYPING.out.hla_result.splitCsv(sep: '\t', header: false, skip:1)
                            .map { row ->
                            def firstColumn = row[0]
                            def otherColumns = row[1..-1].collect { value ->
                                def hla = value.replaceAll('P', '')
                                return "HLA-${hla}"
                            }
                            return [firstColumn, otherColumns]
                            }


        ch_annot_vcf.map { meta,vcf ->
                    def id = meta.id
                    return [id,vcf]
                    }
                    .set{ ch_vcf }

        pvacseq_geno.combine(ch_vcf, by:0)
                    .map{ meta, hla, vcf ->
                    def name = [:]
                    name.id = meta.toString()
                    name.caller = "mutect2"
                    [name, hla, vcf]}
                    .set{ch_hla_vcf}

        EPITOPE_PREDICTION(
            QUANTIFY_SALMON.out.results,
            params.input,
            params.design,
            params.batch,
            ch_hla_vcf,
            params.callers,
            params.neoantigen_epitope1_lengths,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fasta_fai.map { [ [:], it ] },
            PREPARE_GENOME.out.gtf
        )
        ch_versions         = ch_versions.mix( EPITOPE_PREDICTION.out.versions)
    }

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
