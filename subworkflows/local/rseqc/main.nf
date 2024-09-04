// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { GTF2BED                         } from '../../../modules/local/gtf2bed'
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
    gtf           // channel: [ gtf ]

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
    // Run RSeQC tin.py
    //
    tin_txt = Channel.empty()

    RSEQC_TIN(bam_bai, bed)
    tin_txt      = RSEQC_TIN.out.txt
    ch_versions  = ch_versions.mix(RSEQC_TIN.out.versions.first())

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()

    RSEQC_READDISTRIBUTION(bam_bai, bed)
    readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
    ch_versions          = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())

    //
    // Run RSeQC geneBody_coverage.py
    //
    genebodycoverage_rscript   = Channel.empty()

    RSEQC_GENEBODYCOVERAGE(bam_bai, bed)
    genebodycoverage_rscript   = RSEQC_GENEBODYCOVERAGE.out.rscript
    ch_versions                = ch_versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions.first())

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_all     = Channel.empty()
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()

    RSEQC_JUNCTIONSATURATION(bam_bai, bed)
    junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
    junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
    junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
    ch_versions                = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())

    //
    // Run RSeQC Tin Summary
    //
    tin_summary = Channel.empty()

    RSEQC_TINSUMMARY(tin_txt)
    tin_summary                = RSEQC_TINSUMMARY.out.summary_txt
    ch_versions                = ch_versions.mix(RSEQC_TINSUMMARY.out.versions.first())

    /*
    // Run RSeQC Gene Body Coverage Plot
    //
    genebodycoverage_plot    = Channel.empty()

    RSEQC_GENEBODYCOVERAGEPLOT(genebodycoverage_rscript)
    genebodycoverage_plot      = RSEQC_GENEBODYCOVERAGEPLOT.out.png_curves
    ch_versions                = ch_versions.mix(RSEQC_GENEBODYCOVERAGEPLOT.out.versions.first())*/

    //
    // Run RSeQC Read Distribution Matrix
    //
    readdistribution_matrix = Channel.empty()

    RSEQC_READDISTRIBUTIONMATRIX(readdistribution_txt.map{ it[1] }.collect())
    readdistribution_matrix    = RSEQC_READDISTRIBUTIONMATRIX.out.matrix_tab
    ch_versions                = ch_versions.mix(RSEQC_READDISTRIBUTIONMATRIX.out.versions.first())

    emit:

    tin_txt                         // channel: [ val(meta), txt ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    genebodycoverage_rscript

    junctionsaturation_all          // channel: [ val(meta), {pdf, r} ]
    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    tin_summary                     // channel: [ txt ]

    readdistribution_matrix             // channel: [ tab ]

    versions = ch_versions          // channel: [ versions.yml ]
}

