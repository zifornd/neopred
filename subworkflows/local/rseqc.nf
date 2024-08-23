// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { RSEQC_TIN                 } from '../modules/nf-core/rseqc/tin/main'
include { RSEQC_READDISTRIBUTION    } from '../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_JUNCTIONSATURATION  } from '../modules/nf-core/rseqc/junctionsaturation/main'
include { RSEQC_GENEBODYCOVERAGE    } from '../modules/local/rseqc_genebodycoverage'

workflow RSEQC {

    take:
    bam_bai       // channel: [ val(meta), [ bam, bai ] ]
    bed           // channel: [ genome.bed ]

    main:

    bam = bam_bai.map{ [ it[0], it[1] ] }

    ch_versions = Channel.empty()

    //
    // Run RSeQC tin.py
    //
    tin_txt = Channel.empty()

    RSEQC_TIN(bam_bai, bed)
    tin_txt      = RSEQC_TIN.out.txt
    versions    = versions.mix(RSEQC_TIN.out.versions.first())

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()

    RSEQC_READDISTRIBUTION(bam, bed)
    readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
    versions             = versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())

    //
    // Run RSeQC geneBody_coverage.py
    //
    genebodycoverage_rscript   = Channel.empty()

    RSEQC_GENEBODYCOVERAGE(bam, bed)
    genebodycoverage_rscript   = RSEQC_GENEBODYCOVERAGE.out.rscript
    versions                   = versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions.first())

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_all     = Channel.empty()
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()

    RSEQC_JUNCTIONSATURATION(bam, bed)
    junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
    junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
    junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
    versions                   = versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())

    emit:

    tin_txt                         // channel: [ val(meta), txt ]
    
    readdistribution_txt            // channel: [ val(meta), txt ]

    genebodycoverage_rscript

    junctionsaturation_all          // channel: [ val(meta), {pdf, r} ]
    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    versions = ch_versions          // channel: [ versions.yml ]
}

