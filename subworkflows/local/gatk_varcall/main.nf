//
// Variant Identification and Filtering with GATK4
//

include { GATK4_MUTECT2                     } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES          } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL   } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_CALCULATECONTAMINATION      } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS           } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_SELECTVARIANTS              } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { BCFTOOLS_VIEW                     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX                    } from '../../../modules/nf-core/bcftools/index/main'
include { GATK4_COUNTVARIANTS               } from '../../../modules/local/gatk4/countvariants'

workflow VARIANT_CALLINGFILTERING {

    take:
    bam
    bai
    fasta                   // channel: /path/to/genome.fa
    fai                     // channel: /path/to/genome.fai
    dict                    // channel: /path/to/genome.dict
    germline_resource       // channel: /path/to/reference.vcf.gz
    germline_resource_tbi   // channel: /path/to/genome.vcf.gz.tbi
    pon                     // channel: /path/to/genome.vcf.gz
    pon_tbi                 // channel: /path/to/genome.vcf.gz.tbi
    dbsnp                   // channel: /path/to/reference.vcf.gz
    dbsnp_tbi               // channel: /path/to/reference.vcf.gz.tbi
    pileup_vcf
    pileup_vcftbi

    main:

    ch_versions = Channel.empty()

    //
    // Variant Identification using Mutect2
    //

    bam_bai_int = bam
        .join(bai, by: [0])
        .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }


    if (params.variant_filtering) {

        GATK4_MUTECT2 (
        bam_bai_int,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        pon,
        pon_tbi
        )

        variants            =  GATK4_MUTECT2.out.vcf
        variants_tbi        =  GATK4_MUTECT2.out.tbi
        variants_stats      =  GATK4_MUTECT2.out.stats
        ch_versions         =  ch_versions.mix(GATK4_MUTECT2.out.versions.first())
        ch_f1r2             =  GATK4_MUTECT2.out.f1r2



        bam_bai_int = bam
            .join(bai, by: [0])
            .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }

        ch_contaminationtable = Channel.empty()
        ch_tumoursegmentationtable = Channel.empty()

        GATK4_GETPILEUPSUMMARIES(
            bam_bai_int,
            fasta,
            fai,
            dict,
            pileup_vcf,
            pileup_vcftbi
        )

        ch_table = GATK4_GETPILEUPSUMMARIES.out.table
            .map{ meta, table -> [ meta, table, [] ] }
        ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions.first())

        GATK4_LEARNREADORIENTATIONMODEL ( ch_f1r2 )
        ch_artifactprior    = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior
        ch_versions         = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())

        GATK4_CALCULATECONTAMINATION ( ch_table )
        ch_contaminationtable       =  GATK4_CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta - meta.subMap('num_intervals'), cont ] }
        ch_tumoursegmentationtable  = GATK4_CALCULATECONTAMINATION.out.segmentation.map{ meta, seg -> [ meta - meta.subMap('num_intervals'), seg ] }
        ch_versions                 = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first())

        vcf_to_filter = variants.join(variants_tbi, failOnDuplicate: true, failOnMismatch: true)
                            .join(variants_stats, failOnDuplicate: true, failOnMismatch: true)
                            .join(ch_artifactprior, failOnDuplicate: true, failOnMismatch: true)
                            .join(ch_tumoursegmentationtable)
                            .join(ch_contaminationtable)
                        .map{ meta, vcf, tbi, stats, orientation, seg, cont -> [ meta, vcf, tbi, stats, orientation, seg, cont, [] ] }

        GATK4_FILTERMUTECTCALLS (
            vcf_to_filter,
            fasta,
            fai,
            dict
        )

        ch_filtered_vcf = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
                        .map{ meta, vcf, tbi -> [ meta, vcf, tbi ]}
        ch_filtered_tbi = GATK4_FILTERMUTECTCALLS.out.tbi
        ch_versions     = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())

        BCFTOOLS_VIEW (
            ch_filtered_vcf,
            [],
            [],
            []
        )
        ch_bcfview_vcf = BCFTOOLS_VIEW.out.vcf
        ch_versions     = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

        BCFTOOLS_INDEX (
            ch_bcfview_vcf
        )
        ch_bcfview_tbi  = BCFTOOLS_INDEX.out.tbi
        ch_versions     = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

        ch_bcfview_index = ch_bcfview_vcf.join(ch_bcfview_tbi, failOnDuplicate: true, failOnMismatch: true)
                        .map{ meta, vcf, tbi -> [ meta, vcf, tbi, [] ]}

        GATK4_SELECTVARIANTS (
            ch_bcfview_index
        )
        ch_selected_vcf = GATK4_SELECTVARIANTS.out.vcf
        ch_selected_tbi = GATK4_SELECTVARIANTS.out.tbi
        ch_versions     = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

        GATK4_COUNTVARIANTS (
            ch_selected_vcf,
            ch_selected_tbi
        )
        ch_variants_counts =  GATK4_COUNTVARIANTS.out.counts
        ch_versions     = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

    } else {

        GATK4_MUTECT2 (
        bam_bai_int,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        pon,
        pon_tbi
        )

        variants            =  GATK4_MUTECT2.out.vcf
        variants_tbi        =  GATK4_MUTECT2.out.tbi
        variants_stats      =  GATK4_MUTECT2.out.stats
        ch_versions         =  ch_versions.mix(GATK4_MUTECT2.out.versions.first())
        ch_f1r2             =  GATK4_MUTECT2.out.f1r2

        ch_variants_vcf = GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
                        .map{ meta, vcf, tbi -> [ meta, vcf, tbi ]}
        ch_filtered_tbi = GATK4_MUTECT2.out.tbi
        ch_versions     = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

        BCFTOOLS_VIEW (
            ch_variants_vcf,
            [],
            [],
            []
        )
        ch_bcfview_vcf = BCFTOOLS_VIEW.out.vcf
        ch_versions     = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

        BCFTOOLS_INDEX (
            ch_bcfview_vcf
        )
        ch_bcfview_tbi  = BCFTOOLS_INDEX.out.tbi
        ch_versions     = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

        ch_bcfview_index = ch_bcfview_vcf.join(ch_bcfview_tbi, failOnDuplicate: true, failOnMismatch: true)
                        .map{ meta, vcf, tbi -> [ meta, vcf, tbi, [] ]}

        GATK4_SELECTVARIANTS (
            ch_bcfview_index
        )
        ch_selected_vcf = GATK4_SELECTVARIANTS.out.vcf
        ch_selected_tbi = GATK4_SELECTVARIANTS.out.tbi
        ch_versions     = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

        GATK4_COUNTVARIANTS (
            ch_selected_vcf,
            ch_selected_tbi
        )
        ch_variants_counts =  GATK4_COUNTVARIANTS.out.counts
        ch_versions     = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        //ch_versions = Channel.empty()
        ch_table                       =   Channel.empty()
        ch_artifactprior               =   Channel.empty()
        ch_contaminationtable          =   Channel.empty()
        ch_tumoursegmentationtable     =   Channel.empty()
        ch_filtered_vcf                =   Channel.empty()
        ch_filtered_tbi                =   Channel.empty()
    }
        emit:
        table                       =   ch_table
        artifactprior               =   ch_artifactprior
        contaminationtable          =   ch_contaminationtable
        tumour_segmentationtable    =   ch_tumoursegmentationtable
        variants                    =   variants
        filtered_vcf                =   ch_filtered_vcf
        filtered_tbi                =   ch_filtered_tbi
        selected_vcf                =   ch_selected_vcf
        selected_tbi                =   ch_selected_tbi
        versions                    =   ch_versions
}
