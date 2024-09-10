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
include { GATK4_COUNTVARIANTS as GATK4_SNP_COUNTS              } from '../../../modules/local/gatk4_countvariants'

workflow VARIANT_CALLINGFILTERING {

    take:
    bam
    bai
    pileup_vcf
    pileup_vcftbi
    f1r2
    variants
    variants_tbi
    variants_stats
    fasta
    fai
    dict

    main:

    ch_versions = Channel.empty()

    //
    // Variant Identification using Mutect2
    //
    GATK4_MUTECT2 (
    bam,
    fasta,
    fai,
    dict,
    germline_resource,
    pon,
    germline_resource_tbi,
    pon_tbi
    )

    ch_variants_rna_vcf      =  GATK4_MUTECT2.out.vcf
    ch_variants_rna_tbi      =  GATK4_MUTECT2.out.tbi
    ch_variants_stats        =  GATK4_MUTECT2.out.stats
    ch_versions              =  ch_versions.mix(GATK4_MUTECT2.out.versions.first())
    ch_f1r2                  =  GATK4_MUTECT2.out.f1r2

    ch_contaminationtable = Channel.empty()
    ch_tumoursegmentationtable = Channel.empty()
    GATK4_GETPILEUPSUMMARIES(
        bam,
        bai,
        pileup_vcf,
        pileup_vcftbi
    )
    ch_table    = GATK4_GETPILEUPSUMMARIES.out.table
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions.first())

    GATK4_LEARNREADORIENTATIONMODEL ( f1r2 )
    ch_artifactprior    = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior
    ch_versions         = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())

    GATK4_CALCULATECONTAMINATION ( ch_table )
    ch_contaminationtable       =  GATK4_CALCULATECONTAMINATION.out.contamination
    ch_tumoursegmentationtable  = GATK4_CALCULATECONTAMINATION.out.segmentation
    ch_versions                 = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first())


    GATK4_FILTERMUTECTCALLS (
        variants,
        variants_tbi,
        variants_stats,
        ch_artifactprior,
        ch_contaminationtable,
        ch_tumoursegmentationtable,
        fasta,
        fai,
        dict
    )
    ch_filtered_vcf = GATK4_FILTERMUTECTCALLS.out.vcf
    ch_filtered_tbi = GATK4_FILTERMUTECTCALLS.out.tbi
    ch_versions     = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())

    BCFTOOLS_VIEW (
        ch_filtered_vcf,
        ch_filtered_tbi,
    )
    ch_bcfview_vcf  = BCFTOOLS_VIEW.out.vcf
    ch_versions     = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    BCFTOOLS_INDEX (
        ch_bcfview_vcf
    )
    ch_bcfview_tbi  = BCFTOOLS_INDEX.out.tbi
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    GATK4_SELECTVARIANTS (
        ch_bcfview_vcf,
        ch_bcfview_tbi,
        fasta,
        fai,
        dict
    )
    ch_selected_vcf = GATK4_SELECTVARIANTS.out.vcf
    ch_selected_tbi = GATK4_SELECTVARIANTS.out.tbi
    ch_versions     = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

    GATK4_SNP_COUNTS ( ch_selected_vcf, ch_selected_tbi )

    emit:
    table                       =   ch_table
    artifactprior               =   ch_artifactprior
    contaminationtable          =   ch_contaminationtable
    tumour_segmentationtable    =   ch_tumoursegmentationtable
    filtered_vcf                =   ch_filtered_vcf
    filtered_tbi                =   ch_filtered_tbi
    selected_vcf                =   ch_selected_vcf
    selected_tbi                =   ch_selected_tbi
    versions                    =   ch_versions
}
