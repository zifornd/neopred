include { GTF2BED                         } from '../../../modules/local/gtf2bed'
include { RSEQC_BAM_DOWNSAMPLING          } from '../../../modules/local/rseqc/bam_downsampling/main'
include { SAMTOOLS_INDEX                  } from '../../../modules/nf-core/samtools/index'
include { BEDTOOLS_INTERSECT              } from '../../../modules/nf-core/bedtools/intersect/main.nf'
include { SAMTOOLS_INDEX as DOWN_HK_INDEX } from '../../../modules/nf-core/samtools/index'
include { RSEQC_TIN                       } from '../../../modules/nf-core/rseqc/tin/main'
include { RSEQC_READDISTRIBUTION          } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_JUNCTIONSATURATION        } from '../../../modules/nf-core/rseqc/junctionsaturation/main'
include { RSEQC_GENEBODYCOVERAGE          } from '../../../modules/local/rseqc/genebodycoverage'
include { RSEQC_TINSUMMARY                } from '../../../modules/local/rseqc/tinsummary'
//include { RSEQC_GENEBODYCOVERAGEPLOT      } from '../../../modules/local/rseqc/genebodycoverageplot'
include { RSEQC_READDISTRIBUTIONMATRIX    } from '../../../modules/local/rseqc/readdistributionmatrix'

workflow RSEQC {

    take:
    bam_bai       // channel: [ val(meta), [ bam, bai ] ]
    samtools_stats
    gtf           // channel: [ gtf ]
    hk_bed

    main:

    bam = bam_bai.map{ [ it[0], it[1] ] }

    ch_versions = Channel.empty()

    //
    //Module: GTF2BED
    //
    GTF2BED (params.gtf)
    bed              = GTF2BED.out.bed
    ch_versions      = ch_versions.mix(GTF2BED.out.versions)


    //
    // Run RSeQC bam downsampling
    //
    down_bam = Channel.empty()

    RSEQC_BAM_DOWNSAMPLING(bam_bai, samtools_stats)
    down_bam     = RSEQC_BAM_DOWNSAMPLING.out.downsample_bam
    ch_versions  = ch_versions.mix(RSEQC_BAM_DOWNSAMPLING.out.versions.first())

    //
    // Module: Samtools index
    //
    SAMTOOLS_INDEX ( RSEQC_BAM_DOWNSAMPLING.out.downsample_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    RSEQC_BAM_DOWNSAMPLING.out.downsample_bam
    .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
    .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
    .map {
        meta, bam, bai, csi ->
            if (bai) {
                [ meta, bam, bai ]
            } else {
                [ meta, bam, csi ]
            }
    }
    .set { ch_down_bam_bai }
    println ch_down_bam_bai

    //
    // Module: Bedtools intersect
    //

    if ( params.hk_bed ) {

        BEDTOOLS_INTERSECT ( ch_down_bam_bai, hk_bed )
        hk_bam_bai  = BEDTOOLS_INTERSECT.out.hk_bam
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

        //
        // Module: Samtools index
        //
        DOWN_HK_INDEX ( BEDTOOLS_INTERSECT.out.hk_bam )
        ch_versions = ch_versions.mix(DOWN_HK_INDEX.out.versions)

        BEDTOOLS_INTERSECT.out.hk_bam
        .join(DOWN_HK_INDEX.out.bai, by: [0], remainder: true)
        .join(DOWN_HK_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_hk_bam_bai }
        println ch_hk_bam_bai
    }

    //
    // Run RSeQC tin.py
    //


    if ( params.hk_bed ) {

        tin_txt = Channel.empty()

        RSEQC_TIN(ch_hk_bam_bai, params.hk_bed)
        tin_txt      = RSEQC_TIN.out.txt
        ch_versions  = ch_versions.mix(RSEQC_TIN.out.versions.first())

    } else {

        tin_txt = Channel.empty()

        RSEQC_TIN(ch_down_bam_bai, bed)
        tin_txt      = RSEQC_TIN.out.txt
        ch_versions  = ch_versions.mix(RSEQC_TIN.out.versions.first())

    }

/*
    tin_txt = Channel.empty()

    RSEQC_TIN(bam_bai, bed)
    tin_txt      = RSEQC_TIN.out.txt
    ch_versions  = ch_versions.mix(RSEQC_TIN.out.versions.first())
*/
    //
    // Run RSeQC read_distribution.py
    //

    if ( params.hk_bed ) {

        readdistribution_txt = Channel.empty()

        RSEQC_READDISTRIBUTION(ch_hk_bam_bai, params.hk_bed)
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        ch_versions          = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())

    }

    else {
        readdistribution_txt = Channel.empty()

        RSEQC_READDISTRIBUTION(ch_down_bam_bai, bed)
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        ch_versions          = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    //
    // Run RSeQC geneBody_coverage.py
    //

    if ( params.hk_bed ) {

        genebodycoverage_rscript   = Channel.empty()

        RSEQC_GENEBODYCOVERAGE(ch_hk_bam_bai, params.hk_bed)
        genebodycoverage_rscript   = RSEQC_GENEBODYCOVERAGE.out.rscript
        ch_versions                = ch_versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions.first())

    }

    else {

        genebodycoverage_rscript   = Channel.empty()

        RSEQC_GENEBODYCOVERAGE(ch_down_bam_bai, bed)
        genebodycoverage_rscript   = RSEQC_GENEBODYCOVERAGE.out.rscript
        ch_versions                = ch_versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions.first())

    }

    //
    // Run RSeQC junction_saturation.py
    //

    if ( params.hk_bed ) {

    junctionsaturation_all     = Channel.empty()
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()

    RSEQC_JUNCTIONSATURATION(ch_hk_bam_bai, params.hk_bed)
    junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
    junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
    junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
    ch_versions                = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())

    }

    else {

    junctionsaturation_all     = Channel.empty()
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()

    RSEQC_JUNCTIONSATURATION(ch_down_bam_bai, bed)
    junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
    junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
    junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
    ch_versions                = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())

    }

    //
    // Run RSeQC Tin Summary
    //
    tin_summary = Channel.empty()

    RSEQC_TINSUMMARY(tin_txt)
    tin_summary                = RSEQC_TINSUMMARY.out.summary_txt
    ch_versions                = ch_versions.mix(RSEQC_TINSUMMARY.out.versions.first())

/*
    //
    // Run RSeQC Gene Body Coverage Plot
    //
    genebodycoverage_plot    = Channel.empty()

    RSEQC_GENEBODYCOVERAGEPLOT(genebodycoverage_rscript)
    genebodycoverage_plot      = RSEQC_GENEBODYCOVERAGEPLOT.out.png_curves
    ch_versions                = ch_versions.mix(RSEQC_GENEBODYCOVERAGEPLOT.out.versions.first())

*/

    //
    // Run RSeQC Read Distribution Matrix
    //
    readdistribution_matrix = Channel.empty()

    RSEQC_READDISTRIBUTIONMATRIX(readdistribution_txt.map{ it[1] }.collect())
    readdistribution_matrix    = RSEQC_READDISTRIBUTIONMATRIX.out.matrix_tab
    ch_versions                = ch_versions.mix(RSEQC_READDISTRIBUTIONMATRIX.out.versions.first())

    emit:
    down_bam_bai               = ch_down_bam_bai
    //hk_bam_bai                 = ch_hk_bam_bai

    tin_txt                         // channel: [ val(meta), txt ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    genebodycoverage_rscript

    junctionsaturation_all          // channel: [ val(meta), {pdf, r} ]
    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    tin_summary                     // channel: [ txt ]

    readdistribution_matrix         // channel: [ tab ]

    versions                   = ch_versions          // channel: [ versions.yml ]
}

